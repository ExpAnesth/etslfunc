# etslfunc

is a collection of functions dealing with 'time stamp lists' (tsl), 'extended time stamp lists' (etsl), and variants thereof. They are formats for storing time information of discrete neuronal events. See doc/primer_time_stamp_lists.rtf. Essentially all graphical user interface projects in the ExpAn organization which analyze electrophysiological data rely on the code in this repository.

Matlab toolboxes required (for some functions): 
* Signal Processing
* Statistics and Machine Learning

## General note on repositories in the ExpAnesth organization
The code in these repositories provides basic tools for the analysis of electrophysiological time series to members of the Section of Experimental Anesthesiology, Department of Anesthesiology, University Hospital of Tuebingen. Except where noted, code was written by Harald Hentschke. It has been designed primarily for in-house use by individuals who were instructed on its scope and limitations. Also, a substantial proportion of the code has been developed and extended over a time span of >10 years. In detail,

* the implementation of algorithms reflects the evolution of Matlab itself, that is, code that had been initially developed on older versions of Matlab does not necessarily feature newer techniques such as the new automatic array expansion as introduced in Matlab Release 2016b
* nonetheless, all code has been tested to run on Matlab R2018b
* while most m-files contain ample comments, documentation exists only for a few repositories
* checks of user input are implemented to varying degrees
* the code will be improved, updated and documented when and where the need arises