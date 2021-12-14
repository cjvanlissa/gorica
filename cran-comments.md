Note that this package was archived from CRAN:
  "X-CRAN-Comment: Archived on 2021-09-23 as check problems were not
    corrected in time."

We apologize; we had fixed the check problems locally, but forgot to upload the
update to CRAN.

# Version 0.1.2

* Fixed gorica.t_test as per Prof. dr. Brian Ripley's notification
* Bugfix to complement
* Add warnings if penalty terms are < 0
* Add warning if user tries to compute complement for more than 1 hypothesis
* Add Leonard Vanbrabant as contributor

## Test environments

* Local Windows 10 (x64 and x32, build 15063), R 4.1.0
* rhub::check_for_cran()
* devtools::release()
* devtools::oldrelease()
* devtools::devel()

## R CMD check results

OK
