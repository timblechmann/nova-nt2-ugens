HyperbolSaturation : UGen {
	*ar { |sig level|
		^this.multiNew( 'audio', sig, level )
	}
	*kr { |sig level|
		^this.multiNew( 'control', sig, level )
	}
}

ParabolSaturation : HyperbolSaturation {
}

LeakDC2 : MultiOutUGen {

	*ar { arg left, right, cutoff;
		^this.multiNew( 'audio', left, right, cutoff )
	}

	*kr { arg left, right, cutoff;
		^this.multiNew( 'control', left, right, cutoff )
	}

	init { arg ... theInputs;
		inputs = theInputs;
		channels = [ OutputProxy(rate, this, 0), OutputProxy(rate, this, 1) ];
		^channels
	}

	checkInputs { ^this.checkNInputs(2) }
}

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

NovaLowPass2_4th    : NovaLowPass2 {}
NovaHighPass2_4th   : NovaLowPass2 {}
NovaBandPass2_4th   : NovaLowPass2 {}
NovaBandReject2_4th : NovaLowPass2 {}
NovaAllPass2_4th    : NovaLowPass2 {}

NovaLowPass2_2  : NovaLowPass2 {
    *ar { arg left, right, cutoffLeft, cutoffRight, qLeft = 0.70710678118655, qRight = 0.70710678118655;
        ^this.multiNew( 'audio', left, right, cutoffLeft, cutoffRight, qLeft, qRight );
    }

    *kr { arg left, right, cutoffLeft, cutoffRight, qLeft = 0.70710678118655, qRight = 0.70710678118655;
        ^this.multiNew( 'control', left, right, cutoffLeft, cutoffRight, qLeft, qRight );
    }
}

NovaLowPass2_2_4th : NovaLowPass2_2 {}

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
