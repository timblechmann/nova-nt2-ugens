NovaLeakDC : PureUGen {
    *ar { arg sig, cutoff;
        ^this.multiNew( 'audio', sig, cutoff )
    }

    *kr { arg sig, cutoff;
        ^this.multiNew( 'control', sig, cutoff )
    }
}

NovaLeakDC2 : PureMultiOutUGen {
    *ar { arg left, right, cutoff;
        ^this.multiNew( 'audio', left, right, cutoff )
    }

    *kr { arg left, right, cutoff;
        ^this.multiNew( 'control', left, right, cutoff )
    }

    init { arg ... theInputs;
        inputs = theInputs;
        channels = (0..1).collect( OutputProxy(rate, this, _) );
        ^channels
    }

    checkInputs { ^this.checkNInputs(2) }
}

NovaLeakDC4 : PureMultiOutUGen {
    *ar { arg s0, s1, s2, s3, cutoff;
        ^this.multiNew( 'audio', s0, s1, s2, s3, cutoff )
    }

    *kr { arg s0, s1, s2, s3, cutoff;
        ^this.multiNew( 'control', s0, s1, s2, s3, cutoff )
    }

    init { arg ... theInputs;
        inputs = theInputs;
        channels = (0..3).collect( OutputProxy(rate, this, _) );
        ^channels
    }

    checkInputs { ^this.checkNInputs(4) }
}

NovaLeakDC8 : PureMultiOutUGen {
    *ar { arg s0, s1, s2, s3, s4, s5, s6, s7, cutoff;
        ^this.multiNew( 'audio', s0, s1, s2, s3, s4, s5, s6, s7, cutoff )
    }

    *kr { arg s0, s1, s2, s3, s4, s5, s6, s7, cutoff;
        ^this.multiNew( 'control', s0, s1, s2, s3, s4, s5, s6, s7, cutoff )
    }

    init { arg ... theInputs;
        inputs = theInputs;
        channels = (0..7).collect( OutputProxy(rate, this, _) );
        ^channels
    }

    checkInputs { ^this.checkNInputs(8) }
}


NovaIntegrator2 : MultiOutUGen {
    *ar { arg left, right, cutoff;
        ^this.multiNew( 'audio', left, right, cutoff )
    }

    *kr { arg left, right, cutoff;
        ^this.multiNew( 'control', left, right, cutoff )
    }

    init { arg ... theInputs;
        inputs = theInputs;
        channels = (0..2).collect( OutputProxy(rate, this, _) );
        ^channels
    }

    checkInputs { ^this.checkNInputs(2) }
}

NovaIntegrator4 : MultiOutUGen {
    *ar { arg s0, s1, s2, s3, cutoff;
        ^this.multiNew( 'audio', s0, s1, s2, s3, cutoff )
    }

    *kr { arg s0, s1, s2, s3, cutoff;
        ^this.multiNew( 'control', s0, s1, s2, s3, cutoff )
    }

    init { arg ... theInputs;
        inputs = theInputs;
        channels = (0..3).collect( OutputProxy(rate, this, _) );
        ^channels
    }

    checkInputs { ^this.checkNInputs(4) }
}

NovaIntegrator8 : MultiOutUGen {
    *ar { arg s0, s1, s2, s3, s4, s5, s6, s7, cutoff;
        ^this.multiNew( 'audio', s0, s1, s2, s3, s4, s5, s6, s7, cutoff )
    }

    *kr { arg s0, s1, s2, s3, s4, s5, s6, s7, cutoff;
        ^this.multiNew( 'control', s0, s1, s2, s3, s4, s5, s6, s7, cutoff )
    }

    init { arg ... theInputs;
        inputs = theInputs;
        channels = (0..7).collect( OutputProxy(rate, this, _) );
        ^channels
    }

    checkInputs { ^this.checkNInputs(8) }
}


//////////////////////

NovaLowPass : PureUGen {
    *ar { arg sig, cutoff = 440, q = 0.70710678118655;
		^this.multiNew( 'audio', sig, cutoff, q )
    }

    *kr { arg sig, cutoff = 440, q = 0.70710678118655;
        ^this.multiNew( 'control', sig, cutoff, q )
    }

    checkInputs { ^this.checkNInputs(1) }
}

NovaHighPass   : NovaLowPass {}
NovaBandPass   : NovaLowPass {}
NovaBandReject : NovaLowPass {}
NovaAllPass    : NovaLowPass {}

// 2-channel
NovaLowPass2 : PureMultiOutUGen {
    *ar { arg left, right, cutoff = 440, q = 0.70710678118655;
        ^this.multiNew( 'audio', left, right, cutoff, q )
    }

    *kr { arg left, right, cutoff = 440, q = 0.70710678118655;
        ^this.multiNew( 'control', left, right, cutoff, q )
    }

    init { arg ... theInputs;
        inputs = theInputs;
        channels = [ OutputProxy(rate, this, 0), OutputProxy(rate, this, 1) ];
        ^channels
    }

    checkInputs { ^this.checkNInputs(2) }
}

NovaHighPass2   : NovaLowPass2 {}
NovaBandPass2   : NovaLowPass2 {}
NovaBandReject2 : NovaLowPass2 {}
NovaAllPass2    : NovaLowPass2 {}


// 2-channel, separate controls
NovaLowPass2_2 : PureMultiOutUGen {
    *ar { arg left, right, cutoffLeft = 440, cutoffRight = 440, qLeft = 0.70710678118655, qRight = 0.70710678118655;
        ^this.multiNew( 'audio', left, right, cutoffLeft, cutoffRight, qLeft, qRight )
    }

    *kr { arg left, right, cutoffLeft = 440, cutoffRight = 440, qLeft = 0.70710678118655, qRight = 0.70710678118655;
        ^this.multiNew( 'control', left, right, cutoffLeft, cutoffRight, qLeft, qRight )
    }

    init { arg ... theInputs;
        inputs = theInputs;
        channels = [ OutputProxy(rate, this, 0), OutputProxy(rate, this, 1) ];
        ^channels
    }

    checkInputs { ^this.checkNInputs(2) }
}


NovaHighPass2_2   : NovaLowPass2_2 {}
NovaBandPass2_2   : NovaLowPass2_2 {}
NovaBandReject2_2 : NovaLowPass2_2 {}
NovaAllPass2_2    : NovaLowPass2_2 {}


// 4-channel

NovaLowPass4 : PureMultiOutUGen {
    *ar { arg s0, s1, s2, s3, cutoff = 440, q = 0.70710678118655;
        ^this.multiNew( 'audio', s0, s1, s2, s3, cutoff, q )
    }

    *kr { arg s0, s1, s2, s3, cutoff = 440, q = 0.70710678118655;
        ^this.multiNew( 'control', s0, s1, s2, s3, cutoff, q )
    }

    init { arg ... theInputs;
        inputs = theInputs;
        channels = [ OutputProxy(rate, this, 0), OutputProxy(rate, this, 1), OutputProxy(rate, this, 2), OutputProxy(rate, this, 3) ];
        ^channels
    }

    checkInputs { ^this.checkNInputs(4) }
}

NovaHighPass4   : NovaLowPass4 {}
NovaBandPass4   : NovaLowPass4 {}
NovaBandReject4 : NovaLowPass4 {}
NovaAllPass4    : NovaLowPass4 {}



// 4-channel

NovaLowPass4_4 : PureMultiOutUGen {
    *ar { arg s0, s1, s2, s3, cutoff0 = 440, cutoff1 = 440, cutoff2 = 440, cutoff3 = 440, q0 = 0.70710678118655, q1 = 0.70710678118655, q2 = 0.70710678118655, q3 = 0.70710678118655;
        ^this.multiNew( 'audio', s0, s1, s2, s3, cutoff0, cutoff1, cutoff2, cutoff3, q0, q1, q2, q3 )
    }

	*kr { arg s0, s1, s2, s3, cutoff0 = 440, cutoff1 = 440, cutoff2 = 440, cutoff3 = 440, q0 = 0.70710678118655, q1 = 0.70710678118655, q2 = 0.70710678118655, q3 = 0.70710678118655;
        ^this.multiNew( 'control', s0, s1, s2, s3, cutoff0, cutoff1, cutoff2, cutoff3, q0, q1, q2, q3 )
    }

    init { arg ... theInputs;
        inputs = theInputs;
        channels = [ OutputProxy(rate, this, 0), OutputProxy(rate, this, 1), OutputProxy(rate, this, 2), OutputProxy(rate, this, 3) ];
        ^channels
    }

    checkInputs { ^this.checkNInputs(4) }
}

NovaHighPass4_4   : NovaLowPass4_4 {}
NovaBandPass4_4   : NovaLowPass4_4 {}
NovaBandReject4_4 : NovaLowPass4_4 {}
NovaAllPass4_4    : NovaLowPass4_4 {}
