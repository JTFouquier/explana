install_r_packages <- function(package_list) {

  package_list <- package_list[!(package_list %in%
  installed.packages()[, "Package"])]
  if (length(package_list) > 0) {
    install.packages(package_list,
    repos = "http://cran.us.r-project.org")
  }

}