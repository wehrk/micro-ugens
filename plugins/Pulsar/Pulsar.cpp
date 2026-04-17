/*
    Pulsar - Pulsar Synthesis UGen for SuperCollider
    
    A pulsar synthesizer with:
    - Polyphonic pulsarets (overlapping when width > 1)
    - Subsample-accurate phase accumulation
    - Polynomial pulsaret shapes: t^a * (1-t)^b
    - Formant control via sine cycles within pulsaret envelope
    - Kuramoto-style phase coupling with multiple inputs
    - Hard/soft sync modes
    - Phase output for coupling multiple oscillators
    
    Shape parameters (shapeA, shapeB):
    - Peak position in envelope = a/(a+b)
    - Peak sharpness = a+b (higher = narrower = brighter spectrum)
    - Asymmetry = a/b ratio (affects attack/decay character)
    
    Spectral relationship:
    - Narrower pulsaret (high a,b) → more harmonics, brighter
    - Wider pulsaret (low a,b) → fewer harmonics, darker
    - Asymmetric (a≠b) → changes phase relationships, affects transient quality
    
    Synchronization modes:
    - Kuramoto (mode 0): continuous phase attraction/repulsion
      Computes sum of sin(2π(θ_j - θ_self)) for all input phases
    - Hard sync (mode 1): reset phase on sync trigger (uses first input only)
    - Soft sync (mode 2): reset only if phase > threshold (uses first input only)
    
    Outputs:
    - OUT(0): audio signal
    - OUT(1): current phase [0-1] for coupling to other oscillators
    
    Input indices:
    0: freq           - pulsar train frequency in Hz
    1: formantRatio   - formant/fundamental ratio (>= 0.125)
    2: width          - duty cycle [0-2], >1 creates overlapping pulsarets
    3: shapeA         - polynomial exponent a [0.1-10]
    4: shapeB         - polynomial exponent b [0.1-10]
    5: mask           - gate: 0 = silence, >0 = output
    6: syncStrength   - coupling strength (unclamped, >1 for chaos)
    7: syncMode       - 0=kuramoto, 1=hard, 2=soft
    8: syncThreshold  - for soft sync, phase threshold [0-1]
    9: phase          - initial phase [0-1] (only read at init)
    10: maxPulsarets  - maximum overlapping pulsarets (1-16, default 4)
    11: numSyncInputs - number of sync phase inputs (1-32, default 1)
    12+: syncPhases   - phase inputs from other oscillators [0-1] each
*/

#include "SC_PlugIn.h"
#include <cmath>
#include <cstring>

static InterfaceTable *ft;

static const double TWOPI = 6.283185307179586;
static const double PI = 3.141592653589793;
static const int MAX_PULSARETS = 16;
static const int MAX_SYNC_INPUTS = 32;
static const int FIRST_SYNC_INPUT = 12;

// ============================================================
// Helper functions
// ============================================================

static inline double pulsar_clamp(double x, double lo, double hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}

static inline float pulsar_clampf(float x, float lo, float hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}

static inline float pulsar_sanitize(float x) {
    return std::isfinite(x) ? x : 0.0f;
}

// ============================================================
// Pulsaret state
// ============================================================

struct Pulsaret {
    bool active;
    double phase;           // envelope phase [0, 1)
    double phaseInc;        // phase increment per sample
    double formantRatio;    // cached formant ratio
    double shapeA;          // cached shape
    double shapeB;
    double invPeakVal;      // 1 / peak of t^a * (1-t)^b (cached at trigger)
};

// ============================================================
// Polynomial pulsaret envelope (peak is pre-normalized by caller)
// ============================================================

static inline double computeEnvelope(double t, double a, double b, double invPeakVal) {
    if (t <= 0.0 || t >= 1.0) return 0.0;
    double raw = std::pow(t, a) * std::pow(1.0 - t, b);
    return raw * invPeakVal;
}

// ============================================================
// Pulsar Unit struct
// ============================================================

struct Pulsar : public Unit {
    double m_phase;           // main train phase [0, 1)
    float m_prevSyncPhase;    // for edge detection in hard/soft sync (first input)
    
    int m_maxPulsarets;       // maximum simultaneous pulsarets
    int m_numSyncInputs;      // number of sync phase inputs
    Pulsaret m_pulsarets[MAX_PULSARETS];
    
    double m_sampleRate;
    double m_sampleDur;
};

extern "C" {
    void Pulsar_Ctor(Pulsar *unit);
    void Pulsar_next(Pulsar *unit, int nSamples);
}

// ============================================================
// Constructor
// ============================================================

void Pulsar_Ctor(Pulsar *unit) {
    unit->m_sampleRate = SAMPLERATE;
    unit->m_sampleDur = SAMPLEDUR;
    
    // Initial phase from input 9
    float initPhase = ZIN0(9);
    initPhase = pulsar_clampf(initPhase, 0.0f, 1.0f);
    unit->m_phase = (double)initPhase;
    unit->m_prevSyncPhase = 0.0f;
    
    // Max pulsarets from input 10
    int maxP = (int)ZIN0(10);
    if (maxP < 1) maxP = 4;
    if (maxP > MAX_PULSARETS) maxP = MAX_PULSARETS;
    unit->m_maxPulsarets = maxP;
    
    // Number of sync inputs from input 11
    int numSync = (int)ZIN0(11);
    if (numSync < 1) numSync = 1;
    if (numSync > MAX_SYNC_INPUTS) numSync = MAX_SYNC_INPUTS;
    unit->m_numSyncInputs = numSync;
    
    // Initialize pulsaret slots
    for (int i = 0; i < MAX_PULSARETS; ++i) {
        unit->m_pulsarets[i].active = false;
        unit->m_pulsarets[i].phase = 0.0;
    }
    
    SETCALC(Pulsar_next);
    Pulsar_next(unit, 1);
}

// ============================================================
// Find free pulsaret slot
// ============================================================

static int findFreePulsaret(Pulsar *unit) {
    // Find inactive slot
    for (int i = 0; i < unit->m_maxPulsarets; ++i) {
        if (!unit->m_pulsarets[i].active) return i;
    }
    // If all active, steal the one closest to finishing
    int best = 0;
    double highestPhase = 0.0;
    for (int i = 0; i < unit->m_maxPulsarets; ++i) {
        if (unit->m_pulsarets[i].phase > highestPhase) {
            highestPhase = unit->m_pulsarets[i].phase;
            best = i;
        }
    }
    return best;
}

// ============================================================
// Trigger new pulsaret
// ============================================================

static void triggerPulsaret(Pulsar *unit, double width, double formantRatio,
                            double shapeA, double shapeB, double freq,
                            double subsampleOffset) {
    int idx = findFreePulsaret(unit);
    Pulsaret *p = &unit->m_pulsarets[idx];

    p->active = true;

    // Pulsaret completes when phase reaches 1
    // Duration = width / freq seconds
    // phaseInc = 1 / (duration * sampleRate) = freq / (width * sampleRate)
    p->phaseInc = freq / (width * unit->m_sampleRate);

    // Subsample-accurate start: pulsaret was triggered subsampleOffset samples
    // before the next processing step. At its first processed sample it should
    // already be subsampleOffset * phaseInc into the envelope.
    double startPhase = subsampleOffset * p->phaseInc;
    if (startPhase < 0.0) startPhase = 0.0;
    if (startPhase > 0.999) startPhase = 0.999;
    p->phase = startPhase;

    // Anti-alias the carrier: pulsaret carrier frequency is
    // freq * formantRatio / width (phase advances at freq/(width*SR) per
    // sample, carrier is sin(2π·phase·formantRatio)). Clamp formantRatio so
    // the carrier stays below 0.49 * SR.
    double carrierHz = freq * formantRatio / width;
    double nyquistLimit = unit->m_sampleRate * 0.49;
    if (carrierHz > nyquistLimit && freq > 1e-9) {
        formantRatio = nyquistLimit * width / freq;
        if (formantRatio < 0.125) formantRatio = 0.125;
    }

    p->formantRatio = formantRatio;
    p->shapeA = shapeA;
    p->shapeB = shapeB;

    // Cache the envelope peak normalization
    double peakT = shapeA / (shapeA + shapeB);
    double peakVal = std::pow(peakT, shapeA) * std::pow(1.0 - peakT, shapeB);
    p->invPeakVal = (peakVal > 1e-12) ? (1.0 / peakVal) : 0.0;
}

// ============================================================
// DSP calculation
// ============================================================

void Pulsar_next(Pulsar *unit, int nSamples) {
    float *outAudio = OUT(0);
    float *outPhase = OUT(1);
    
    // Input pointers for fixed inputs
    float *freqIn = IN(0);
    float *formantRatioIn = IN(1);
    float *widthIn = IN(2);
    float *shapeAIn = IN(3);
    float *shapeBIn = IN(4);
    float *maskIn = IN(5);
    float *syncStrengthIn = IN(6);
    float *syncModeIn = IN(7);
    float *syncThresholdIn = IN(8);
    
    // Check rates for fixed inputs
    int freqRate = INRATE(0);
    int formantRate = INRATE(1);
    int widthRate = INRATE(2);
    int shapeARate = INRATE(3);
    int shapeBRate = INRATE(4);
    int maskRate = INRATE(5);
    int syncStrengthRate = INRATE(6);
    int syncModeRate = INRATE(7);
    int syncThresholdRate = INRATE(8);
    
    double phase = unit->m_phase;
    float prevSyncPhase = unit->m_prevSyncPhase;
    double sampleDur = unit->m_sampleDur;
    double sampleRate = unit->m_sampleRate;
    int numSyncInputs = unit->m_numSyncInputs;
    
    for (int i = 0; i < nSamples; ++i) {
        // Get parameters
        float freq = (freqRate == calc_FullRate) ? freqIn[i] : freqIn[0];
        float formantRatio = (formantRate == calc_FullRate) ? formantRatioIn[i] : formantRatioIn[0];
        float width = (widthRate == calc_FullRate) ? widthIn[i] : widthIn[0];
        float shapeA = (shapeARate == calc_FullRate) ? shapeAIn[i] : shapeAIn[0];
        float shapeB = (shapeBRate == calc_FullRate) ? shapeBIn[i] : shapeBIn[0];
        float mask = (maskRate == calc_FullRate) ? maskIn[i] : maskIn[0];
        float syncStrength = (syncStrengthRate == calc_FullRate) ? syncStrengthIn[i] : syncStrengthIn[0];
        float syncMode = (syncModeRate == calc_FullRate) ? syncModeIn[i] : syncModeIn[0];
        float syncThreshold = (syncThresholdRate == calc_FullRate) ? syncThresholdIn[i] : syncThresholdIn[0];
        
        // Sanitize and clamp
        freq = pulsar_sanitize(freq);
        if (freq < 0.001f) freq = 0.001f;
        if (freq > sampleRate * 0.5f) freq = (float)(sampleRate * 0.5);
        
        formantRatio = pulsar_sanitize(formantRatio);
        if (formantRatio < 0.125f) formantRatio = 0.125f;
        if (formantRatio > 64.0f) formantRatio = 64.0f;
        
        width = pulsar_clampf(pulsar_sanitize(width), 0.01f, 2.0f);
        shapeA = pulsar_clampf(pulsar_sanitize(shapeA), 0.1f, 10.0f);
        shapeB = pulsar_clampf(pulsar_sanitize(shapeB), 0.1f, 10.0f);
        mask = pulsar_sanitize(mask);
        syncStrength = pulsar_sanitize(syncStrength);
        syncThreshold = pulsar_clampf(pulsar_sanitize(syncThreshold), 0.0f, 1.0f);
        int syncModeInt = (int)pulsar_clampf(pulsar_sanitize(syncMode), 0.0f, 2.0f);
        
        // Base phase increment
        double phaseInc = freq * sampleDur;
        
        // Get first sync phase (used for hard/soft sync edge detection)
        float *firstSyncIn = IN(FIRST_SYNC_INPUT);
        int firstSyncRate = INRATE(FIRST_SYNC_INPUT);
        float firstSyncPhase = (firstSyncRate == calc_FullRate) ? firstSyncIn[i] : firstSyncIn[0];
        firstSyncPhase = pulsar_sanitize(firstSyncPhase);
        
        // Apply synchronization
        if (syncModeInt == 0 && std::abs(syncStrength) > 0.0001f) {
            // Kuramoto-style coupling: sum pairwise coupling terms
            // dθ/dt = ω + (K/N) * Σ sin(θ_j - θ)
            
            double couplingSum = 0.0;
            
            for (int s = 0; s < numSyncInputs; ++s) {
                float *syncIn = IN(FIRST_SYNC_INPUT + s);
                int syncRate = INRATE(FIRST_SYNC_INPUT + s);
                float syncPhase = (syncRate == calc_FullRate) ? syncIn[i] : syncIn[0];
                syncPhase = pulsar_sanitize(syncPhase);
                
                // Phase difference wrapped to [-0.5, 0.5]
                double phaseDiff = (double)syncPhase - phase;
                if (phaseDiff > 0.5) phaseDiff -= 1.0;
                if (phaseDiff < -0.5) phaseDiff += 1.0;
                
                // Basic Kuramoto coupling term
                couplingSum += std::sin(TWOPI * phaseDiff);
            }
            
            // Normalize by number of inputs
            double coupling = couplingSum / (double)numSyncInputs;
            
            double K = (double)syncStrength;
            double absK = std::abs(K);
            
            // Add nonlinear terms for chaos (based on total coupling strength)
            if (absK > 0.5) {
                // Second harmonic - use average phase difference for this
                double avgPhaseDiff = 0.0;
                for (int s = 0; s < numSyncInputs; ++s) {
                    float *syncIn = IN(FIRST_SYNC_INPUT + s);
                    int syncRate = INRATE(FIRST_SYNC_INPUT + s);
                    float syncPhase = (syncRate == calc_FullRate) ? syncIn[i] : syncIn[0];
                    double pd = (double)pulsar_sanitize(syncPhase) - phase;
                    if (pd > 0.5) pd -= 1.0;
                    if (pd < -0.5) pd += 1.0;
                    avgPhaseDiff += pd;
                }
                avgPhaseDiff /= (double)numSyncInputs;
                
                coupling += 0.5 * std::sin(2.0 * TWOPI * avgPhaseDiff);
                
                if (absK > 1.0) {
                    coupling += 0.3 * std::sin(3.0 * TWOPI * avgPhaseDiff);
                    double cubicTerm = avgPhaseDiff * avgPhaseDiff * avgPhaseDiff * 4.0;
                    coupling += cubicTerm * (absK - 1.0);
                }
                if (absK > 1.5) {
                    double selfMod = std::sin(TWOPI * phase * 2.0) * (absK - 1.5) * 0.5;
                    coupling += selfMod;
                }
            }
            
            double couplingInc = K * coupling * freq * sampleDur;
            double maxInc = phaseInc * 4.0;
            couplingInc = std::tanh(couplingInc / maxInc) * maxInc;
            
            phaseInc += couplingInc;
        }
        else if (syncModeInt == 1 || syncModeInt == 2) {
            // Hard/soft sync - uses first sync input only
            if (prevSyncPhase > 0.8f && firstSyncPhase < 0.2f) {
                // Master wrapped between samples i-1 and i. Linearly interpolate
                // across the wrap (dSync is the wrap-corrected phase delta) to
                // find how far into the current sample the wrap occurred.
                double dSync = (double)firstSyncPhase + 1.0 - (double)prevSyncPhase;
                double subOff = 0.5;
                if (dSync > 1e-9) {
                    double fracToWrap = (1.0 - (double)prevSyncPhase) / dSync;
                    subOff = 1.0 - fracToWrap;
                    if (subOff < 0.0) subOff = 0.0;
                    if (subOff > 1.0) subOff = 1.0;
                }

                bool doSync = (syncModeInt == 1) ||
                              (syncModeInt == 2 && phase > (double)syncThreshold);

                if (doSync) {
                    // Own phase should equal subOff * phaseInc at this sample
                    // (as if reset subOff samples ago).
                    phase = subOff * phaseInc;
                    if (mask >= 0.5f) {
                        triggerPulsaret(unit, width, formantRatio, shapeA, shapeB,
                                        freq, subOff);
                    }
                }
            }
        }
        
        prevSyncPhase = firstSyncPhase;
        
        // Sum all active pulsarets
        double sample = 0.0;
        
        for (int p = 0; p < unit->m_maxPulsarets; ++p) {
            Pulsaret *pulsaret = &unit->m_pulsarets[p];
            if (!pulsaret->active) continue;

            // Compute envelope (peak normalization cached at trigger)
            double env = computeEnvelope(pulsaret->phase, pulsaret->shapeA,
                                         pulsaret->shapeB, pulsaret->invPeakVal);

            // Compute carrier
            double carrier = std::sin(TWOPI * pulsaret->phase * pulsaret->formantRatio);

            sample += env * carrier;

            // Advance pulsaret phase
            pulsaret->phase += pulsaret->phaseInc;

            if (pulsaret->phase >= 1.0) {
                pulsaret->active = false;
            }
        }

        outAudio[i] = (float)sample;
        outPhase[i] = (float)phase;

        // Advance main phase
        phase += phaseInc;

        // Check for new pulsaret trigger
        if (phase >= 1.0) {
            // Subsample-accurate trigger: overshoot tells us how far past the
            // wrap we've advanced, so the new pulsaret was triggered
            // (overshoot / phaseInc) samples ago.
            double overshoot = phase - std::floor(phase);
            double subOff = 0.0;
            if (phaseInc > 1e-12) {
                subOff = (phase - 1.0) / phaseInc;
                if (subOff < 0.0) subOff = 0.0;
                if (subOff > 1.0) subOff = 1.0;
            }
            phase = overshoot;

            // Trigger new pulsaret if not masked
            if (mask >= 0.5f) {
                triggerPulsaret(unit, width, formantRatio, shapeA, shapeB,
                                freq, subOff);
            }
        }
        else if (phase < 0.0) {
            phase += 1.0 - std::floor(phase);
        }
    }
    
    unit->m_phase = phase;
    unit->m_prevSyncPhase = prevSyncPhase;
}

// ============================================================
// Plugin interface
// ============================================================

PluginLoad(Pulsar) {
    ft = inTable;
    DefineSimpleUnit(Pulsar);
}
