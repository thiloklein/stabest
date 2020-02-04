# stabest
Estimate preferences in a two-sided school market with stability restrictions

*IMPORTANT*: since the function has the same name as the package, the automatic documentation system will mess up the .Rd docs and reference both stabit and stabit-package to the same files. So, the workflow to create the documentation is:

1. devtools::document()
2. I manually need to delete \alias{stabit-package} in ./man/stabit.Rd, and \alias{stabit-package} in ./man/stabit-package.Rd.
3. build the package w/o enabling Roxygen functionality
