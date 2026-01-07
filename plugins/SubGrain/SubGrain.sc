SubGrain : MultiOutUGen {
    // Subsample-accurate granular synthesis
    // 
    // Features:
    // - Subsample-accurate grain scheduling (smooth frequency sweeps)
    // - Internal delay line for live input granulation
    // - External buffer for sample-based granulation (MONO ONLY)
    // - Density/regularity morphing (periodic to stochastic)
    // - Anti-aliased interpolation (cubic, 6-point, sinc)
    // - Morphable envelope (percussive to hanning to reversed)
    // - Multichannel panning with spread
    //
    // All per-grain parameters are sampled at grain birth.
    
    *ar {
        arg numChannels = 2,        // number of output channels (scalar)
            maxGrains = 32,         // maximum simultaneous grains (scalar)
            maxDelayTime = 1.0,     // delay line length in seconds (scalar)
            
            input = 0,              // audio input for delay line (ar)
            buffer = -1,            // buffer number, MONO ONLY (kr)
            
            bufferMix = 1.0,        // 0 = delay only, 1 = buffer only (kr/ar)
            
            density = 10,           // grains per second (kr/ar)
            regularity = 1.0,       // 0 = stochastic, 1 = periodic (kr/ar)
            
            duration = 0.1,         // grain duration in seconds (kr/ar)
            rate = 1.0,             // playback rate, per grain (kr/ar)
            
            position = 0.5,         // position [0-1] in buffer/delay (kr/ar)
                                        // For delay: 0=oldest(max delay), 1=newest(min delay)
                                        // Use ~0.9-0.99 for responsive live granulation
            positionSpread = 0,     // random position spread (kr/ar)
            
            envShape = 0.5,         // 0=percussive, 0.5=hanning, 1=reversed (kr/ar)
            
            pan = 0.5,              // pan position [0-1] (kr/ar)
            panSpread = 0,          // random pan spread (kr/ar)
            
            interpolation = 0;      // 0=cubic, 1=6-point, 2=sinc (kr)
        
        ^this.multiNew('audio',
            numChannels, maxGrains, maxDelayTime,
            input, buffer,
            bufferMix,
            density, regularity,
            duration, rate,
            position, positionSpread,
            envShape,
            pan, panSpread,
            interpolation
        )
    }
    
    init { arg ... theInputs;
        inputs = theInputs;
        ^this.initOutputs(theInputs[0].asInteger, rate);
    }
    
    checkInputs {
        // numChannels must be scalar
        if (inputs[0].rate != 'scalar') {
            ^"numChannels must be scalar"
        };
        // maxGrains must be scalar
        if (inputs[1].rate != 'scalar') {
            ^"maxGrains must be scalar"
        };
        // maxDelayTime must be scalar
        if (inputs[2].rate != 'scalar') {
            ^"maxDelayTime must be scalar"
        };
        ^this.checkValidInputs
    }
}
