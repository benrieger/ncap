#! /bin/bash

# Copyright (C) 2010, Christian Meesters (meesters@uni-mainz.de)
# This code is part of the ncap distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# cleanup the doc directory after a LaTeX run

\rm -f ./doc/*.blg
\rm -f ./doc/*.out
\rm -f ./doc/*.aux
\rm -f ./doc/*~
\rm -f ./doc/*.log
\rm -f ./doc/*.toc
\rm -f ./doc/*.bbl
\rm -f ./doc/*.backup

