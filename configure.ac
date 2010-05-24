AC_INIT("DESCRIPTION")

AC_CHECK_SIZEOF([unsigned long])
AC_CHECK_LIB([z], [gzeof],,[AC_MSG_ERROR([need libz with gzeof])])
## AC_CHECK_LIB([z], [gzopen],,[AC_MSG_ERROR([need libz with gzopen])])
## AC_CHECK_LIB([z], [gzclose],,[AC_MSG_ERROR([need libz with gzclose])])
## AC_CHECK_LIB([z], [gzgets],,[AC_MSG_ERROR([need libz with gzgets])])
## AC_CHECK_LIB([z], [gzgetc],,[AC_MSG_ERROR([need libz with gzgetc])])

AC_OUTPUT(src/Makevars)
