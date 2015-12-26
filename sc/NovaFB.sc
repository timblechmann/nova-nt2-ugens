NovaFBIn : MultiOutUGen {
	*ar {|numberOfChannels = 1|
		^this.multiNew('audio', numberOfChannels);
	}

	init { arg argNumChannels ... theInputs;
		inputs = [argNumChannels];
		^this.initOutputs(argNumChannels + 1, rate)
	}
}

NovaFBOut : UGen {
	*ar {|inUgen, numberOfChannels, signal|
		var args = [inUgen, numberOfChannels] ++ signal;

		^this.multiNew('audio', *args)
	}
}

NovaFBNode {
	var <channels, fbInNode, fbOutNode;

	*new {|channelCount|
		^super.newCopyArgs(channelCount)
	}

	read {
		var ret = NovaFBIn.ar(channels);

		if (fbInNode.notNil) {
			Error("NovaFBNode: only one read call is allowed").throw;
		};

		fbInNode = ret[0];
		^ret[1..]
	}

	write {|args|
		if (args.size != channels) {
			Error("NovaFBNode: channel count mismatch % %".format(args.size, channels)).throw;
		};

		if (fbOutNode.notNil) {
			Error("NovaFBNode: write called multiple times").throw;
		};

		fbOutNode = NovaFBOut.ar(fbInNode, channels, args);
	}
}
