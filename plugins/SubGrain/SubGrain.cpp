/*
    SubGrain - Subsample-Accurate Granular Synthesis UGen
    
    Features:
    - Subsample-accurate grain triggering for smooth frequency control
    - Density/regularity triggering (periodic to stochastic)
    - Internal delay line for live input + external buffer (MONO ONLY)
    - Anti-aliased playback (cubic/6-point/sinc interpolation)
    - Morphable envelope (percussive to hanning to reversed)
    - Multichannel equal-power panning
    - All parameters sampled per grain (including rate)
    
    IMPORTANT: Buffer must be mono. Multi-channel buffers will be rejected.
*/

#include "SC_PlugIn.h"
#include <cmath>
#include <cstring>

static InterfaceTable *ft;

static const double TWOPI = 6.283185307179586;
static const double PI = 3.141592653589793;
static const int SINC_TAPS = 8;

// ============================================================
// Sinc lookup table (normalized sinc: sin(π·w) / (π·w))
// ============================================================
// Range [0, SINC_LUT_MAX_ARG], sinc is even so negative args fold.
// Max arg = 32 so halfTaps ≤ 32 with cutoff ≤ 1 is representable.
static const int SINC_LUT_MAX_ARG = 32;
static const int SINC_LUT_RES = 1024;
static const int SINC_LUT_SIZE = SINC_LUT_MAX_ARG * SINC_LUT_RES + 1;
static float g_sincLut[SINC_LUT_SIZE];
static bool g_sincLutInit = false;

static void initSincLut() {
    if (g_sincLutInit) return;
    g_sincLut[0] = 1.0f;
    for (int i = 1; i < SINC_LUT_SIZE; ++i) {
        double x = (double)i / (double)SINC_LUT_RES;
        g_sincLut[i] = (float)(std::sin(PI * x) / (PI * x));
    }
    g_sincLutInit = true;
}

static inline float sincLut(double w) {
    double aw = std::abs(w);
    if (aw >= (double)SINC_LUT_MAX_ARG) return 0.0f;
    double idx = aw * (double)SINC_LUT_RES;
    int i0 = (int)idx;
    if (i0 >= SINC_LUT_SIZE - 1) return 0.0f;
    double frac = idx - (double)i0;
    float a = g_sincLut[i0];
    float b = g_sincLut[i0 + 1];
    return a + (float)frac * (b - a);
}

// Grain state
struct Grain {
    bool active;
    
    double phase;           // envelope phase [0, 1]
    double phaseInc;        // phase increment per sample
    
    double readPos;         // samples read since grain start
    double readInc;         // rate (samples per sample)
    float bufferMix;        // blend between delay and buffer
    
    double bufferStart;     // starting position in buffer (samples)
    double delayStart;      // starting delay time (samples)

    float envShape;         // 0=percussive, 0.5=hanning, 1=reversed

    // Source buffer NUMBER captured at grain birth, so a buffer-number swap
    // mid-stream only affects grains born after it — live grains finish reading
    // the buffer they started on (no mid-grain content jump / click on a
    // double-buffer swap). We store the number (not the data pointer) and resolve
    // the live SndBuf each read, so an in-place re-alloc of the same buffer can't
    // leave a grain reading freed memory. -1 = no buffer.
    int bufnum;

    float* panAmps;         // pre-computed per-channel amplitudes
};

struct SubGrain : public Unit {
    int m_numChannels;
    int m_maxGrains;
    float m_maxDelayTime;
    
    Grain* m_grains;
    float* m_panAmpStorage;
    
    float* m_delayLine;
    int m_delayFrames;
    int m_delayWritePos;
    
    double m_triggerPhase;
    
    uint32 m_randState;
    
    double m_sampleRate;
    double m_sampleDur;
    
    // Buffer cache - refreshed each block
    const float* m_bufData;
    int m_bufFrames;
    int m_bufnum;            // current source buffer number (-1 = none/invalid); captured per grain at birth

    float* m_outputAccum;
    float* m_inputCache;     // block-size scratch for the audio input
    int m_maxBlockSize;

    bool m_allocFailed;
    int m_lastWarnedBufnum;  // avoid spamming the multi-channel warning
};

// Forward declarations
extern "C" {
    void SubGrain_Ctor(SubGrain *unit);
    void SubGrain_Dtor(SubGrain *unit);
    void SubGrain_next(SubGrain *unit, int nSamples);
    void SubGrain_next_failed(SubGrain *unit, int nSamples);
}

// ============================================================
// Helper functions
// ============================================================

static inline double sg_clamp(double x, double lo, double hi) {
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

static inline float sg_clampf(float x, float lo, float hi) { 
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

static inline bool sg_isfinite(float x) {
    return std::isfinite(x);
}

static inline bool sg_isfinite_d(double x) {
    return std::isfinite(x);
}

static inline float sg_sanitize(float x) {
    return sg_isfinite(x) ? x : 0.0f;
}

// Fast wrap using fmod (avoids expensive while loops)
static inline double sg_wrap(double x, double hi) {
    if (hi <= 0.0) return 0.0;
    x = std::fmod(x, hi);
    if (x < 0.0) x += hi;
    return x;
}

static inline float getInput(SubGrain* unit, int index, int sampleIndex) {
    Wire* wire = unit->mInput[index];
    if (wire->mCalcRate == calc_FullRate) {
        return sg_sanitize(wire->mBuffer[sampleIndex]);
    } else {
        return sg_sanitize(wire->mBuffer[0]);
    }
}

static inline double nextRandom(SubGrain* unit) {
    uint32 x = unit->m_randState;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    unit->m_randState = x;
    return (double)x / 4294967296.0;
}

// ============================================================
// Panning
// ============================================================

static void computePanAmps(SubGrain* unit, float pan, float* amps) {
    int numCh = unit->m_numChannels;
    
    if (numCh == 1) {
        amps[0] = 1.0f;
        return;
    }
    
    pan = sg_clampf(pan, 0.0f, 1.0f);
    
    if (numCh == 2) {
        float angle = pan * (float)(PI * 0.5);
        amps[0] = std::cos(angle);
        amps[1] = std::sin(angle);
        return;
    }
    
    // Multichannel: circular panning
    float angle = pan * (float)TWOPI;
    float channelSpacing = (float)TWOPI / numCh;
    
    for (int ch = 0; ch < numCh; ++ch) {
        float channelAngle = ch * channelSpacing;
        float diff = angle - channelAngle;
        
        while (diff > PI) diff -= (float)TWOPI;
        while (diff < -PI) diff += (float)TWOPI;
        
        if (std::abs(diff) < channelSpacing) {
            amps[ch] = std::cos(diff * (float)PI / (2.0f * channelSpacing));
        } else {
            amps[ch] = 0.0f;
        }
    }
}

// ============================================================
// Interpolation
// ============================================================

// 4-point cubic (Hermite)
static float cubicInterp(const float* data, int frames, double pos) {
    if (!data || frames <= 0) return 0.0f;
    if (frames == 1) return sg_sanitize(data[0]);
    if (!sg_isfinite_d(pos)) return 0.0f;
    
    pos = sg_wrap(pos, (double)frames);
    
    int i1 = (int)pos;
    double frac = pos - i1;
    
    int i0 = (i1 - 1 + frames) % frames;
    int i2 = (i1 + 1) % frames;
    int i3 = (i1 + 2) % frames;
    
    float y0 = data[i0], y1 = data[i1], y2 = data[i2], y3 = data[i3];
    
    float c0 = y1;
    float c1 = 0.5f * (y2 - y0);
    float c2 = y0 - 2.5f * y1 + 2.0f * y2 - 0.5f * y3;
    float c3 = 0.5f * (y3 - y0) + 1.5f * (y1 - y2);
    
    float result = ((c3 * (float)frac + c2) * (float)frac + c1) * (float)frac + c0;
    return sg_sanitize(result);
}

// 6-point polynomial (good balance)
static float interp6(const float* data, int frames, double pos) {
    if (!data || frames <= 0) return 0.0f;
    if (frames < 6) return cubicInterp(data, frames, pos);
    if (!sg_isfinite_d(pos)) return 0.0f;
    
    pos = sg_wrap(pos, (double)frames);
    
    int i = (int)pos;
    float x = (float)(pos - i);
    
    int i0 = (i - 2 + frames) % frames;
    int i1 = (i - 1 + frames) % frames;
    int i2 = i % frames;
    int i3 = (i + 1) % frames;
    int i4 = (i + 2) % frames;
    int i5 = (i + 3) % frames;
    
    float y0 = data[i0], y1 = data[i1], y2 = data[i2];
    float y3 = data[i3], y4 = data[i4], y5 = data[i5];
    
    float c0 = y2;
    float c1 = (y3 - y1) * 0.5f;
    float c2 = y1 - y2 * 2.5f + y3 * 2.0f - y4 * 0.5f;
    float c3 = (y4 - y0) * 0.5f + (y1 - y3) * 1.5f;
    float c4 = y0 * 0.25f - y1 + y2 * 1.5f - y3 + y4 * 0.25f;
    float c5 = (y5 - y0) / 12.0f + (y1 - y4) * 0.5f + (y3 - y2) * (5.0f/6.0f);
    
    float result = c0 + x * (c1 + x * (c2 + x * (c3 + x * (c4 + x * c5))));
    return sg_isfinite(result) ? result : 0.0f;
}

// Sinc with anti-aliasing. Uses the shared sinc LUT for both the band-limit
// kernel and the Lanczos window, so no std::sin calls in the inner loop.
// cutoff = min(1, 1/|rate|) for pitch-up.
static float sincInterp(const float* data, int frames, double pos, double rate) {
    if (!data || frames <= 0) return 0.0f;
    if (frames == 1) return sg_sanitize(data[0]);
    if (!sg_isfinite_d(pos) || !sg_isfinite_d(rate)) return 0.0f;

    pos = sg_wrap(pos, (double)frames);

    int center = (int)pos;
    double frac = pos - (double)center;

    double absRate = std::abs(rate);
    if (absRate < 0.001) absRate = 0.001;
    double cutoff = (absRate > 1.0) ? (1.0 / absRate) : 1.0;

    // More taps for narrower filters. Cap matches SINC_LUT_MAX_ARG.
    int effectiveTaps = SINC_TAPS;
    if (cutoff < 0.5) {
        effectiveTaps = (int)(SINC_TAPS / cutoff);
        if (effectiveTaps > 64) effectiveTaps = 64;
    }
    int halfTaps = effectiveTaps / 2;
    if (halfTaps < 2) halfTaps = 2;

    double invHalfTaps = 1.0 / (double)halfTaps;
    float cutoffF = (float)cutoff;

    float sum = 0.0f;
    float norm = 0.0f;

    // Symmetric tap range: t ∈ [-halfTaps+1, halfTaps]
    // so x = t - frac ∈ (-halfTaps, halfTaps], keeping the Lanczos window
    // symmetric and zero at both edges.
    for (int t = -halfTaps + 1; t <= halfTaps; ++t) {
        int idx = ((center + t) % frames + frames) % frames;
        double x = (double)t - frac;

        // Band-limit: cutoff * sinc_norm(cutoff * x)
        float coef = cutoffF * sincLut(x * cutoff);
        // Lanczos window: sinc_norm(x / halfTaps)
        coef *= sincLut(x * invHalfTaps);

        norm += coef;
        sum += data[idx] * coef;
    }

    if (std::abs(norm) > 0.001f) {
        sum /= norm;
    }

    return sg_sanitize(sum);
}

// ============================================================
// Envelope: single parameter morphing
// ============================================================

static float computeEnvelope(double phase, float shape) {
    // shape: 0 = sharp attack + exp decay (percussive)
    //        0.5 = hanning (symmetric)
    //        1 = exp attack + sharp decay (reversed)
    
    if (phase < 0.0 || phase > 1.0 || !sg_isfinite_d(phase)) return 0.0f;
    shape = sg_clampf(shape, 0.0f, 1.0f);
    
    float tilt = 0.02f + shape * 0.96f;
    float decayLen = 1.0f - tilt;
    if (decayLen < 0.02f) decayLen = 0.02f;
    
    float expAmount = std::abs(shape - 0.5f) * 2.0f;
    float expRate = 4.0f + expAmount * 8.0f;
    float expNorm = 1.0f - std::exp(-expRate);
    if (expNorm < 0.001f) expNorm = 0.001f;
    
    float env;
    
    if (phase < tilt) {
        float t = (float)(phase / tilt);
        float cosEnv = 0.5f * (1.0f - std::cos((float)PI * t));
        
        if (shape <= 0.5f) {
            env = cosEnv;
        } else {
            float expEnv = (1.0f - std::exp(-t * expRate)) / expNorm;
            env = cosEnv * (1.0f - expAmount) + expEnv * expAmount;
        }
    } else {
        float t = (float)((phase - tilt) / decayLen);
        float cosEnv = 0.5f * (1.0f + std::cos((float)PI * t));
        
        if (shape >= 0.5f) {
            env = cosEnv;
        } else {
            float expEnv = std::exp(-t * expRate);
            env = cosEnv * (1.0f - expAmount) + expEnv * expAmount;
        }
    }
    
    return sg_clampf(env, 0.0f, 1.0f);
}

// ============================================================
// Grain processing
// ============================================================

static float processGrain(SubGrain* unit, Grain* grain, 
                          const float* delayLine,
                          const float* bufferData, int bufferFrames,
                          int interp) {
    
    float bufferSample = 0.0f;
    float delaySample = 0.0f;
    
    double readOffset = grain->readPos * grain->readInc;
    double absRate = std::abs(grain->readInc);
    
    // Read from buffer
    if (bufferData && bufferFrames > 0) {
        double bufPos = grain->bufferStart + readOffset;
        bufPos = sg_wrap(bufPos, (double)bufferFrames);
        
        if (interp == 2) {
            bufferSample = sincInterp(bufferData, bufferFrames, bufPos, absRate);
        } else if (interp == 1) {
            bufferSample = interp6(bufferData, bufferFrames, bufPos);
        } else {
            bufferSample = cubicInterp(bufferData, bufferFrames, bufPos);
        }
    }
    
    // Read from delay line - delayStart is the START INDEX, advance by readOffset
    if (delayLine && unit->m_delayFrames > 0) {
        double delayPos = grain->delayStart + readOffset;
        delayPos = sg_wrap(delayPos, (double)unit->m_delayFrames);
        
        if (interp == 2) {
            delaySample = sincInterp(delayLine, unit->m_delayFrames, delayPos, absRate);
        } else if (interp == 1) {
            delaySample = interp6(delayLine, unit->m_delayFrames, delayPos);
        } else {
            delaySample = cubicInterp(delayLine, unit->m_delayFrames, delayPos);
        }
    }
    
    float mix = sg_clampf(grain->bufferMix, 0.0f, 1.0f);
    float sample = delaySample * (1.0f - mix) + bufferSample * mix;
    
    grain->readPos += 1.0;
    
    float env = computeEnvelope(grain->phase, grain->envShape);
    sample *= env;
    
    grain->phase += grain->phaseInc;
    
    if (grain->phase >= 1.0) {
        grain->active = false;
    }
    
    return sg_sanitize(sample);
}

// ============================================================
// Grain triggering
// ============================================================

static int findFreeGrain(SubGrain* unit) {
    int bestIdx = 0;
    double highestPhase = -1.0;
    
    for (int i = 0; i < unit->m_maxGrains; ++i) {
        if (!unit->m_grains[i].active) return i;
        if (unit->m_grains[i].phase > highestPhase) {
            highestPhase = unit->m_grains[i].phase;
            bestIdx = i;
        }
    }
    return bestIdx;
}

static void triggerGrain(SubGrain* unit, double subsampleOffset, int sampleIndex, int currentWritePos) {
    if (!unit->m_grains || !unit->m_panAmpStorage) return;
    
    int idx = findFreeGrain(unit);
    Grain* g = &unit->m_grains[idx];
    
    g->active = true;
    g->bufferMix = sg_clampf(getInput(unit, 5, sampleIndex), 0.0f, 1.0f);
    
    float duration = getInput(unit, 8, sampleIndex);
    if (duration <= 0.0f) duration = 0.1f;
    duration = sg_clampf(duration, 0.001f, 10.0f);
    
    float rate = getInput(unit, 9, sampleIndex);
    if (!sg_isfinite(rate)) rate = 1.0f;
    
    g->phaseInc = 1.0 / (duration * unit->m_sampleRate);
    if (g->phaseInc <= 0.0 || !sg_isfinite_d(g->phaseInc)) g->phaseInc = 0.001;
    g->phase = sg_clamp(subsampleOffset * g->phaseInc, 0.0, 0.999);
    
    float position = getInput(unit, 10, sampleIndex);
    float posSpread = getInput(unit, 11, sampleIndex);
    if (!sg_isfinite(position)) position = 0.5f;
    if (!sg_isfinite(posSpread)) posSpread = 0.0f;
    position += (float)(nextRandom(unit) * 2.0 - 1.0) * posSpread;
    position = sg_clampf(position, 0.0f, 1.0f);
    
    // Capture the source buffer NUMBER at grain birth (see Grain::bufnum). The
    // grain reads THIS buffer for its whole life, even if the buffer-number input
    // is swapped to a different buffer while the grain is still playing.
    g->bufnum = unit->m_bufnum;

    // Convert position to buffer samples
    if (unit->m_bufFrames > 0) {
        g->bufferStart = position * (unit->m_bufFrames - 1);
    } else {
        g->bufferStart = 0.0;
    }
    
    // Convert position to delay line START INDEX (not delay time!)
    // position 1 = newest (small delay), position 0 = oldest (max delay)
    //
    // The grain read head advances at `rate` samples/sample while the write
    // head advances at 1 sample/sample. The gap between write and read is
    //     gap(k) = delaySamples + k * (1 - rate)
    // For rate > 1 the gap shrinks, so minDelay must be large enough that the
    // grain never catches the write head. For rate < 0 the gap grows, so
    // maxDelay must leave room before it wraps past the write head.
    double rateD = (double)rate;
    double absRate = std::abs(rateD);
    double durationSamples = (double)duration * unit->m_sampleRate;
    const double interpSafety = 34.0; // max sinc halfTaps + slack
    const double blockSize = (double)unit->mWorld->mFullRate.mBufLength;

    double minDelay = blockSize + interpSafety;
    if (rateD > 1.0) {
        minDelay += durationSamples * (rateD - 1.0);
    }

    double maxDelay = (double)unit->m_delayFrames - interpSafety;
    if (rateD < 0.0) {
        maxDelay -= durationSamples * (1.0 + absRate);
    }

    if (maxDelay < minDelay) {
        // Rate × duration too large for this delay line: grain can't avoid
        // the write head. Clamp to max delay (safer than refusing to play).
        maxDelay = minDelay;
    }

    double delaySamples = (1.0 - (double)position) *
                          ((double)unit->m_delayFrames - interpSafety);
    if (delaySamples < minDelay) delaySamples = minDelay;
    if (delaySamples > maxDelay) delaySamples = maxDelay;
    
    // Store the actual buffer index where this grain starts reading
    // currentWritePos is where we JUST wrote, so currentWritePos-1 is the newest sample
    double startIndex = (double)currentWritePos - delaySamples;
    g->delayStart = sg_wrap(startIndex, (double)unit->m_delayFrames);
    
    // Rate
    g->readInc = (double)rate;
    g->readPos = -subsampleOffset;
    
    // Envelope shape
    g->envShape = getInput(unit, 12, sampleIndex);
    if (!sg_isfinite(g->envShape)) g->envShape = 0.5f;
    g->envShape = sg_clampf(g->envShape, 0.0f, 1.0f);
    
    // Pan
    float pan = getInput(unit, 13, sampleIndex);
    float panSpread = getInput(unit, 14, sampleIndex);
    if (!sg_isfinite(pan)) pan = 0.5f;
    if (!sg_isfinite(panSpread)) panSpread = 0.0f;
    pan += (float)(nextRandom(unit) * 2.0 - 1.0) * panSpread;
    pan = sg_clampf(pan, 0.0f, 1.0f);
    computePanAmps(unit, pan, g->panAmps);
}

// ============================================================
// Constructor
// ============================================================

void SubGrain_Ctor(SubGrain *unit) {
    unit->m_allocFailed = false;
    
    unit->m_numChannels = (int)IN0(0);
    unit->m_maxGrains = (int)IN0(1);
    unit->m_maxDelayTime = IN0(2);
    
    if (unit->m_numChannels < 1) unit->m_numChannels = 1;
    if (unit->m_numChannels > 32) unit->m_numChannels = 32;
    if (unit->m_maxGrains < 1) unit->m_maxGrains = 1;
    if (unit->m_maxGrains > 1024) unit->m_maxGrains = 1024;
    if (unit->m_maxDelayTime < 0.01f) unit->m_maxDelayTime = 0.01f;
    if (unit->m_maxDelayTime > 60.0f) unit->m_maxDelayTime = 60.0f;
    
    unit->m_sampleRate = SAMPLERATE;
    unit->m_sampleDur = 1.0 / unit->m_sampleRate;
    
    // Initialize all pointers to null
    unit->m_delayLine = nullptr;
    unit->m_grains = nullptr;
    unit->m_panAmpStorage = nullptr;
    unit->m_outputAccum = nullptr;
    unit->m_inputCache = nullptr;
    unit->m_maxBlockSize = 0;
    
    // Allocate delay line
    unit->m_delayFrames = (int)(unit->m_maxDelayTime * unit->m_sampleRate) + 64;
    unit->m_delayLine = (float*)RTAlloc(unit->mWorld, unit->m_delayFrames * sizeof(float));
    if (!unit->m_delayLine) {
        unit->m_allocFailed = true;
    } else {
        memset(unit->m_delayLine, 0, unit->m_delayFrames * sizeof(float));
    }
    unit->m_delayWritePos = 0;
    
    // Allocate grain pool
    unit->m_grains = (Grain*)RTAlloc(unit->mWorld, unit->m_maxGrains * sizeof(Grain));
    unit->m_panAmpStorage = (float*)RTAlloc(unit->mWorld, unit->m_maxGrains * unit->m_numChannels * sizeof(float));
    
    if (!unit->m_grains || !unit->m_panAmpStorage) {
        unit->m_allocFailed = true;
    } else {
        memset(unit->m_grains, 0, unit->m_maxGrains * sizeof(Grain));
        memset(unit->m_panAmpStorage, 0, unit->m_maxGrains * unit->m_numChannels * sizeof(float));
        
        for (int i = 0; i < unit->m_maxGrains; ++i) {
            Grain* g = &unit->m_grains[i];
            g->active = false;
            g->phase = 0.0;
            g->phaseInc = 0.001;
            g->readPos = 0.0;
            g->readInc = 1.0;
            g->bufferStart = 0.0;
            g->delayStart = 0.0;
            g->bufferMix = 1.0f;
            g->envShape = 0.5f;
            g->bufnum = -1;
            g->panAmps = unit->m_panAmpStorage + (i * unit->m_numChannels);
            for (int ch = 0; ch < unit->m_numChannels; ++ch) {
                g->panAmps[ch] = (ch == 0) ? 1.0f : 0.0f;
            }
        }
    }
    
    unit->m_triggerPhase = 0.5;  // Avoid immediate trigger
    
    unit->m_randState = (uint32)(unit->mParent->mRGen->irand(0x7FFFFFFF)) | 1;
    for (int i = 0; i < 10; ++i) nextRandom(unit);

    unit->m_bufData = nullptr;
    unit->m_bufFrames = 0;
    unit->m_lastWarnedBufnum = -1;

    int blockSize = unit->mWorld->mFullRate.mBufLength;
    unit->m_maxBlockSize = blockSize;

    // Allocate output accumulator
    unit->m_outputAccum = (float*)RTAlloc(unit->mWorld, blockSize * unit->m_numChannels * sizeof(float));
    if (!unit->m_outputAccum) {
        unit->m_allocFailed = true;
    } else {
        memset(unit->m_outputAccum, 0, blockSize * unit->m_numChannels * sizeof(float));
    }

    // Allocate input cache sized to the actual block length
    unit->m_inputCache = (float*)RTAlloc(unit->mWorld, blockSize * sizeof(float));
    if (!unit->m_inputCache) {
        unit->m_allocFailed = true;
    } else {
        memset(unit->m_inputCache, 0, blockSize * sizeof(float));
    }
    
    // Set calc function based on allocation success
    if (unit->m_allocFailed) {
        SETCALC(SubGrain_next_failed);
        ClearUnitOutputs(unit, 1);
    } else {
        SETCALC(SubGrain_next);
        ClearUnitOutputs(unit, 1);
    }
}

// ============================================================
// Destructor
// ============================================================

void SubGrain_Dtor(SubGrain *unit) {
    if (unit->m_delayLine) RTFree(unit->mWorld, unit->m_delayLine);
    if (unit->m_grains) RTFree(unit->mWorld, unit->m_grains);
    if (unit->m_panAmpStorage) RTFree(unit->mWorld, unit->m_panAmpStorage);
    if (unit->m_outputAccum) RTFree(unit->mWorld, unit->m_outputAccum);
    if (unit->m_inputCache) RTFree(unit->mWorld, unit->m_inputCache);
}

// ============================================================
// Failed allocation - output silence
// ============================================================

void SubGrain_next_failed(SubGrain *unit, int nSamples) {
    ClearUnitOutputs(unit, nSamples);
}

// ============================================================
// Main DSP loop
// ============================================================
// Per-block buffer resolution with supernova reader locks.
//
// scsynth runs buffer commands on the DSP thread, so raw mSndBufs reads are
// safe there. supernova frees/re-fills buffers on a separate non-realtime
// thread WHILE the (multiple) DSP threads run, so every SndBuf must be read
// under its shared lock, held for as long as the data pointer is used. Each
// buffer touched during a block is locked once on first use and released at
// the end of SubGrain_next; the ACQUIRE/RELEASE macros compile to no-ops on
// scsynth, so this path is shared.
// ============================================================

static const int kMaxLockedBufs = 8;

struct SGBufCache {
    int num;
    int bufnums[kMaxLockedBufs];
    SndBuf* bufs[kMaxLockedBufs];
    const float* data[kMaxLockedBufs];
    int frames[kMaxLockedBufs];
};

// Look up bufnum, locking the buffer shared on first touch this block.
// Returns false (data null) for invalid, non-mono or vanished buffers —
// the affected grains just read silence for this block.
static bool sgResolveBuf(SubGrain* unit, SGBufCache* cache, int bufnum,
                         const float** outData, int* outFrames, bool warnNonMono) {
    *outData = nullptr;
    *outFrames = 0;
    if (bufnum < 0 || (uint32)bufnum >= unit->mWorld->mNumSndBufs) return false;
    for (int k = 0; k < cache->num; ++k) {
        if (cache->bufnums[k] == bufnum) {
            *outData = cache->data[k];
            *outFrames = cache->frames[k];
            return cache->data[k] != nullptr;
        }
    }
    if (cache->num >= kMaxLockedBufs) return false;
    SndBuf* buf = unit->mWorld->mSndBufs + bufnum;
    ACQUIRE_SNDBUF_SHARED(buf);
    int k = cache->num++;
    cache->bufnums[k] = bufnum;
    cache->bufs[k] = buf;
    cache->data[k] = nullptr;
    cache->frames[k] = 0;
    if (buf->data && buf->frames > 0) {
        if (buf->channels == 1) {
            cache->data[k] = buf->data;
            cache->frames[k] = buf->frames;
        } else if (warnNonMono && bufnum != unit->m_lastWarnedBufnum) {
            Print("SubGrain: buffer %d has %d channels; only mono "
                  "buffers are supported, ignoring.\n", bufnum, buf->channels);
            unit->m_lastWarnedBufnum = bufnum;
        }
    }
    *outData = cache->data[k];
    *outFrames = cache->frames[k];
    return cache->data[k] != nullptr;
}

static void sgReleaseBufs(SGBufCache* cache) {
    for (int k = 0; k < cache->num; ++k) RELEASE_SNDBUF_SHARED(cache->bufs[k]);
    cache->num = 0;
}

void SubGrain_next(SubGrain *unit, int nSamples) {
    // The input cache, output accumulator, and delay line are all sized to the
    // server's block length, which is constant for the life of the unit.
    if (nSamples > unit->m_maxBlockSize) nSamples = unit->m_maxBlockSize;

    // STEP 1: Cache all inputs BEFORE touching outputs

    // Cache audio input (index 3)
    bool inputIsAudio = (unit->mInput[3]->mCalcRate == calc_FullRate);
    const float* inputBufRaw = unit->mInput[3]->mBuffer;
    float* inputCache = unit->m_inputCache;
    if (inputIsAudio) {
        for (int i = 0; i < nSamples; ++i) {
            inputCache[i] = sg_sanitize(inputBufRaw[i]);
        }
    } else {
        float val = sg_sanitize(inputBufRaw[0]);
        for (int i = 0; i < nSamples; ++i) {
            inputCache[i] = val;
        }
    }
    const float* input = inputCache;
    
    // Refresh buffer cache EVERY BLOCK (buffer can be reallocated). All SndBuf
    // reads this block go through bufCache, which holds each touched buffer's
    // shared lock until the end of this function (supernova; no-op on scsynth).
    SGBufCache bufCache;
    bufCache.num = 0;

    unit->m_bufData = nullptr;
    unit->m_bufFrames = 0;
    unit->m_bufnum = -1;

    float bufnum = sg_sanitize(unit->mInput[4]->mBuffer[0]);
    if (bufnum >= 0.0f) {
        const float* curData;
        int curFrames;
        if (sgResolveBuf(unit, &bufCache, (int)bufnum, &curData, &curFrames, true)) {
            unit->m_bufData = curData;
            unit->m_bufFrames = curFrames;
            unit->m_bufnum = (int)bufnum;
        }
    }
    
    // STEP 2: Clear output accumulator
    if (unit->m_outputAccum) {
        memset(unit->m_outputAccum, 0, nSamples * unit->m_numChannels * sizeof(float));
    }
    
    // STEP 3: Process samples
    int interp = (int)getInput(unit, 15, 0);
    
    for (int i = 0; i < nSamples; ++i) {
        // Write input to delay line
        if (unit->m_delayLine) {
            unit->m_delayLine[unit->m_delayWritePos] = input[i];
            unit->m_delayWritePos = (unit->m_delayWritePos + 1) % unit->m_delayFrames;
        }
        
        // Get trigger parameters
        float density = getInput(unit, 6, i);
        float regularity = getInput(unit, 7, i);
        
        // density <= 0 means no new grains (silence when existing grains finish)
        if (density > 0.0f) {
            regularity = sg_clampf(regularity, 0.0f, 1.0f);
            double phaseInc = density * unit->m_sampleDur;
            
            unit->m_triggerPhase += phaseInc;
            
            // Handle MULTIPLE triggers per sample at high density
            while (unit->m_triggerPhase >= 1.0) {
                double overshoot = unit->m_triggerPhase - 1.0;
                double triggerSubsample = 1.0 - (overshoot / phaseInc);
                triggerSubsample = sg_clamp(triggerSubsample, 0.0, 0.999);
                
                bool doTrigger;
                if (regularity >= 0.999f) {
                    doTrigger = true;
                } else if (regularity <= 0.001f) {
                    doTrigger = false;
                } else {
                    doTrigger = (nextRandom(unit) < regularity);
                }
                
                if (doTrigger) {
                    triggerGrain(unit, triggerSubsample, i, unit->m_delayWritePos);
                }
                
                unit->m_triggerPhase -= 1.0;
            }
            
            // Stochastic triggers (additional random grains when regularity < 1)
            if (regularity < 0.999f) {
                double stochasticProb = (1.0 - regularity) * density * unit->m_sampleDur;
                if (stochasticProb > 0.5) stochasticProb = 0.5;
                if (nextRandom(unit) < stochasticProb) {
                    double subsampleOffset = nextRandom(unit);
                    triggerGrain(unit, subsampleOffset, i, unit->m_delayWritePos);
                }
            }
        }
        
        // Process active grains
        for (int g = 0; g < unit->m_maxGrains; ++g) {
            Grain* grain = &unit->m_grains[g];
            if (!grain->active) continue;
            
            // resolve the grain's captured source buffer through bufCache: safe even
            // if that buffer was re-alloc'd (in-place legacy load) or freed — we read
            // the current valid data under its shared lock, never a dangling pointer.
            const float* gBufData;
            int gBufFrames;
            sgResolveBuf(unit, &bufCache, grain->bufnum, &gBufData, &gBufFrames, false);
            float grainSample = processGrain(unit, grain, unit->m_delayLine,
                                             gBufData, gBufFrames,
                                             interp);
            
            if (unit->m_outputAccum && grain->panAmps) {
                for (int ch = 0; ch < unit->m_numChannels; ++ch) {
                    unit->m_outputAccum[i * unit->m_numChannels + ch] += grainSample * grain->panAmps[ch];
                }
            }
        }
    }
    
    // Done with buffer data — release the shared locks taken this block.
    sgReleaseBufs(&bufCache);

    // STEP 4: Copy accumulator to outputs
    for (int i = 0; i < nSamples; ++i) {
        for (int ch = 0; ch < unit->m_numChannels; ++ch) {
            float val = unit->m_outputAccum ? unit->m_outputAccum[i * unit->m_numChannels + ch] : 0.0f;
            OUT(ch)[i] = sg_sanitize(val);
        }
    }
}

// ============================================================
// Plugin interface
// ============================================================

PluginLoad(SubGrain) {
    ft = inTable;

    initSincLut();

    (*ft->fDefineUnit)("SubGrain",
                       sizeof(SubGrain),
                       (UnitCtorFunc)SubGrain_Ctor,
                       (UnitDtorFunc)SubGrain_Dtor,
                       0);
}
