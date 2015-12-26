////////////////

NovaPanB2D : MultiOutUGen {

    *ar { arg sig, azimut = 0;
        ^this.multiNew( 'audio', sig, azimut )
    }

    *kr { arg sig, azimut;
        ^this.multiNew( 'audio', sig, azimut )
    }

    init { arg ... theInputs;
        inputs = theInputs;
        channels = [ OutputProxy(rate, this, 0), OutputProxy(rate, this, 1), OutputProxy(rate, this, 2)];
        ^channels
    }

    checkInputs { ^this.checkNInputs(1) }
}


DiodeLadderFilter : Filter {
    *ar { arg sig, freq = 440, q = 0.2, feedbackHPF = 1000;
        ^this.multiNew('audio', sig, freq, q, feedbackHPF)
    }
}


DiodeLadderFilter2 : PureMultiOutUGen {
    *ar { arg sig0, sig1, freq = 440, q = 0.2, feedbackHPF = 1000;
        ^this.multiNew('audio', sig0, sig1, freq, q, feedbackHPF)
    }

    init { arg ... theInputs;
        inputs = theInputs;
        channels = [ OutputProxy(rate, this, 0), OutputProxy(rate, this, 1) ];
        ^channels
    }

    checkInputs { ^this.checkNInputs(2) }
}


DiodeLadderFilter4 : PureMultiOutUGen {
    *ar { arg sig0, sig1, sig2, sig3, freq = 440, q = 0.2, feedbackHPF = 1000;
        ^this.multiNew('audio', sig0, sig1, sig2, sig3, freq, q, feedbackHPF)
    }

    init { arg ... theInputs;
        inputs = theInputs;
        channels = [ OutputProxy(rate, this, 0), OutputProxy(rate, this, 1), OutputProxy(rate, this, 2), OutputProxy(rate, this, 3) ];
        ^channels
    }

    checkInputs { ^this.checkNInputs(4) }
}

DiodeLadderFilter4_4 : PureMultiOutUGen {
    *ar { arg sig0, sig1, sig2, sig3, freq0 = 440, freq1 = 440, freq2 = 440, freq3 = 440,
        q0 = 0.2, q1 = 0.2, q2 = 0.2, q3 = 0.2,
        feedbackHPF0 = 1000, feedbackHPF1 = 1000, feedbackHPF2 = 1000, feedbackHPF3 = 1000;
        ^this.multiNew('audio', sig0, sig1, sig2, sig3, freq0, freq1, freq2, freq3,
            q0, q1, q2, q3, feedbackHPF0, feedbackHPF1, feedbackHPF2, feedbackHPF3)
    }

    init { arg ... theInputs;
        inputs = theInputs;
        channels = [ OutputProxy(rate, this, 0), OutputProxy(rate, this, 1), OutputProxy(rate, this, 2), OutputProxy(rate, this, 3) ];
        ^channels
    }

    checkInputs { ^this.checkNInputs(4) }
}
