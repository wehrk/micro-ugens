Pulsar : MultiOutUGen {
    // Pulsar synthesis with polynomial pulsaret shapes and phase coupling
    //
    // Outputs:
    //   [0] - audio signal
    //   [1] - current phase [0-1] for coupling to other Pulsars
    //
    // The pulsaret envelope is t^a * (1-t)^b normalized to peak at 1.0
    //
    // Sync modes:
    //   0 - Kuramoto: sum of sin(θ_j - θ_self) for all inputs
    //   1 - Hard sync: reset on first input's phase wrap
    //   2 - Soft sync: reset only if phase > threshold
    
    *ar {
        arg freq = 100,            // pulsar train frequency (Hz)
            formantRatio = 1,      // formant/fundamental ratio
            width = 0.5,           // duty cycle [0-2], >1 overlaps
            shapeA = 2,            // attack shape exponent [0.1-10]
            shapeB = 2,            // decay shape exponent [0.1-10]
            mask = 1,              // gate: 0=silence, >0=output
            syncPhases = #[0],     // phase input(s) - single value or array
            syncStrength = 0,      // coupling strength (>1 for chaos)
            syncMode = 0,          // 0=kuramoto, 1=hard, 2=soft
            syncThreshold = 0.5,   // threshold for soft sync
            phase = 0,             // initial phase [0-1]
            maxPulsarets = 4,      // max overlapping pulsarets (1-16)
            mul = 1,
            add = 0;
        
        var syncArray = syncPhases.asArray;
        var numSyncInputs = syncArray.size;
        
        ^this.multiNew('audio', freq, formantRatio, width, shapeA, shapeB, 
            mask, syncStrength, syncMode, syncThreshold, phase, maxPulsarets, 
            numSyncInputs, *syncArray)
            .madd(mul, add)
    }
    
    init { arg ... theInputs;
        inputs = theInputs;
        ^this.initOutputs(2, rate);
    }
    
    checkInputs {
        ^this.checkValidInputs
    }
}
