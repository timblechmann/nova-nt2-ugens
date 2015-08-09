HyperbolSaturation : PureUGen {
	*ar { |sig level|
		^this.multiNew( 'audio', sig, level )
	}
	*kr { |sig level|
		^this.multiNew( 'control', sig, level )
	}
}

ParabolSaturation : HyperbolSaturation {
}

PowSaturation : HyperbolSaturation {
}

NovaLeakDC1 : PureUGen {
	*ar { arg sig, cutoff;
		^this.multiNew( 'audio', sig, cutoff )
	}

	*kr { arg sig, cutoff;
		^this.multiNew( 'control', sig, cutoff )
	}

	checkInputs { ^this.checkNInputs(1) }
}

NovaLeakDC2 : MultiOutUGen {
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

NovaLeakDC4 : MultiOutUGen {
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

NovaLeakDC8 : MultiOutUGen {
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

NovaFeedbackAM2 : MultiOutUGen {
        *ar { arg left, right, fb;
                ^this.multiNew( 'audio', left, right, fb )
        }

        *kr { arg left, right, fb;
                ^this.multiNew( 'control', left, right, fb )
        }

        init { arg ... theInputs;
                inputs = theInputs;
                channels = (0..2).collect( OutputProxy(rate, this, _) );
                ^channels
        }

        checkInputs { ^this.checkNInputs(2) }
}

NovaFeedbackAM2_2 : MultiOutUGen {
        *ar { arg left, right, fb0, fb1;
                ^this.multiNew( 'audio', left, right, fb0, fb1 )
        }

        *kr { arg left, right, fb0, fb1;
                ^this.multiNew( 'control', left, right, fb0, fb1 )
        }

        init { arg ... theInputs;
                inputs = theInputs;
                channels = (0..2).collect( OutputProxy(rate, this, _) );
                ^channels
        }

        checkInputs { ^this.checkNInputs(2) }
}

NovaFeedbackAM4_4 : MultiOutUGen {
        *ar { arg s0, s1, s2, s3, fb0, fb1, fb2, fb3;
                ^this.multiNew( 'audio', s0, s1, s2, s3, fb0, fb1, fb2, fb3 )
        }

        *kr { arg s0, s1, s2, s3, fb0, fb1, fb2, fb3;
                ^this.multiNew( 'control', s0, s1, s2, s3, fb0, fb1, fb2, fb3 )
        }

        init { arg ... theInputs;
                inputs = theInputs;
                channels = (0..3).collect( OutputProxy(rate, this, _) );
                ^channels
        }

        checkInputs { ^this.checkNInputs(4) }
}

NovaFeedbackAM8_8 : MultiOutUGen {
        *ar { arg s0, s1, s2, s3, s4, s5, s6, s7, fb0, fb1, fb2, fb3, fb4, fb5, fb6, fb7;
                ^this.multiNew( 'audio', s0, s1, s2, s3, s4, s5, s6, s7, fb0, fb1, fb2, fb3, fb4, fb5, fb6, fb7 )
        }

        *kr { arg s0, s1, s2, s3, s4, s5, s6, s7, fb0, fb1, fb2, fb3, fb4, fb5, fb6, fb7;
                ^this.multiNew( 'control', s0, s1, s2, s3, s4, s5, s6, s7, fb0, fb1, fb2, fb3, fb4, fb5, fb6, fb7 )
        }

        init { arg ... theInputs;
                inputs = theInputs;
                channels = (0..7).collect( OutputProxy(rate, this, _) );
                ^channels
        }

        checkInputs { ^this.checkNInputs(8) }
}


//////////////////////

NovaLowPass : UGen {

	*ar { arg sig, cutoff, q = 0.70710678118655;
		^this.multiNew( 'audio', sig, cutoff, q )
	}

	*kr { arg sig, cutoff, q = 0.70710678118655;
		^this.multiNew( 'control', sig, cutoff, q )
	}

	checkInputs { ^this.checkNInputs(1) }
}

NovaHighPass   : NovaLowPass {}
NovaBandPass   : NovaLowPass {}
NovaBandReject : NovaLowPass {}
NovaAllPass    : NovaLowPass {}



NovaLowPass2 : MultiOutUGen {
	*ar { arg left, right, cutoff, q = 0.70710678118655;
		^this.multiNew( 'audio', left, right, cutoff, q )
	}

	*kr { arg left, right, cutoff, q = 0.70710678118655;
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

NovaLowPass4 : MultiOutUGen {
	*ar { arg s0, s1, s2, s3, cutoff, q = 0.70710678118655;
		^this.multiNew( 'audio', s0, s1, s2, s3, cutoff, q )
	}

	*kr { arg s0, s1, s2, s3, cutoff, q = 0.70710678118655;
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








// NovaLowPass2_2  : NovaLowPass2 {
// 	*ar { arg left, right, cutoffLeft, cutoffRight, qLeft = 0.70710678118655, qRight = 0.70710678118655;
// 		^this.multiNew( 'audio', left, right, cutoffLeft, cutoffRight, qLeft, qRight );
// 	}
//
// 	*kr { arg left, right, cutoffLeft, cutoffRight, qLeft = 0.70710678118655, qRight = 0.70710678118655;
// 		^this.multiNew( 'control', left, right, cutoffLeft, cutoffRight, qLeft, qRight );
// 	}
// }
//
// NovaLowPass2_2_4th : NovaLowPass2_2 {}

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

DiodeLadderFilter4 : Filter {
    *ar { arg sig0, sig1, sig2, sig3, freq = 440, q = 0.2, feedbackHPF = 1000;
        ^this.multiNew('audio', sig0, sig1, sig2, sig3, freq, q, feedbackHPF)
    }
}

DiodeLadderFilter4_4 : Filter {
    *ar { arg sig0, sig1, sig2, sig3, freq0 = 440, freq1 = 440, freq2 = 440, freq3 = 440,
        q0 = 0.2, q1 = 0.2, q2 = 0.2, q3 = 0.2,
        feedbackHPF0 = 1000, feedbackHPF1 = 1000, feedbackHPF2 = 1000, feedbackHPF3 = 1000;
        ^this.multiNew('audio', sig0, sig1, sig2, sig3, freq0, freq1, freq2, freq3,
            q0, q1, q2, q3, feedbackHPF0, feedbackHPF1, feedbackHPF2, feedbackHPF3)
    }
}
