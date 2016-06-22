# transientclumps

----------------

# About #

Data processing for the JCMT Transient Large Program

-----------------

# Running the code #

Add modules to PYTHONPATH or place them in working directory.
Import the modules into your current session or script.

>> from TCGaussclumpsFunctions import *
>> from TCOffsetFunctions import *
>> from TCPrepFunctions import *

First run the preparation function, see the comments in
TCPrepFunctions.py for information on keywords. Images
should be in .sdf format, and run one at a time.

>> prepare_image( *args )

Next run Gaussclumps, see comments in TCGaussclumpsFunctions.py
for more information on keywords. Image used should be output from
prepare_images(), should be in .sdf format and should be run one
at a time.

>> run_gaussclumps( *args )

Next run the matching program, see comments in TCOffsetFunctions.py
for more information on keywords. Catalogs used should be outputs
from run_gaussclumps( *args ), and should be run one at a time
(one catalog run against one reference catalog).

>> source_match( *args )

The output from this function will be the normalization constants 
to align (in position & flux) the target catalog to the reference
catalog. For more information on outputs see comments in 
TCOffsetFunctions.py
