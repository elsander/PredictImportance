# What is this?
This is an R package containing all of the code necessary to replicate
the analyses of Wootton et al. (in review). There are functions that
may also be useful to you even if you aren't trying to replicate our
work. For example, there are functions to generate random networks
using the Cascade Model, Niche Model, and Minimum Potential Niche
Model.

# How do I install it?
This package has a lot of code that is highly specific to a single
publication, so it isn't available on CRAN. Fortunately, it's very
easy to install a package from github.

This installation can all be done from the R command line. First,
install and load the `devtools` package:

    install.packages("devtools")
    library("devtools")

Then simply use the `install_github` function to install this
package. Other R dependencies will automatically be installed as
well. This may take a while, depending on how long it takes for R to
connect to github, and how many dependencies need to be installed.

    devtools::install_github("esander91/PredictImportance")

# How do I find the documentation?
Once you install the package, you can access documentation for any of
the functions using `?function_name`, as with any other package. If
you want a pdf of the documentation for all functions, you can find it
in at `inst/PredictImportanceDocumentation.pdf` in this
repository. This will be available on your computer wherever the
computer stores R packages. Alternatively, you can navigate to the
file on the github website. Click the file and right click the *Raw*
button above the file preview. If you then click *Save Link As...*,
you can save a local copy of the file to your computer.
