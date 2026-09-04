Every file this project changed, exactly as upstream ships it.

They are here so tests/dsk_build_equivalence.c can build an unmodified binary
and compare against it without a network round trip to the upstream
repository. Nothing links against them; build.sh --upstream copies them over
their modified counterparts in a scratch tree and builds that.

Upstream: https://github.com/CoMoS-SA/Reissl_2025.git
Commit:   611ff9cb44348baa55be1bc315eefe2c117ccd44
