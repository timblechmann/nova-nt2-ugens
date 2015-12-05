NovaXFade : PureUGen {
    *ar { arg left, right, pan = 0, amp = 1;
        ^this.multiNew( 'audio', left, right, pan, amp )
    }

    *kr { arg left, right, pan = 0, amp = 1;
        ^this.multiNew( 'control', left, right, pan, amp )
    }
}

LinXFade : PureUGen {
    *ar { arg left, right, pan = 0, amp = 1;
        ^this.multiNew( 'audio', left, right, pan, amp )
    }

    *kr { arg left, right, pan = 0, amp = 1;
        ^this.multiNew( 'control', left, right, pan, amp )
    }
}
