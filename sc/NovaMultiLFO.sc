NovaMultiLFO : UGen {
    *kr { arg mode = 3, freq = 1, phase = 0, min = -1, max = 1;
        ^this.multiNew( 'control', mode, freq, phase, min, max )
    }
}
