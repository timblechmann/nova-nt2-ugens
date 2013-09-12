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
