TAB=`echo -n -e '\t'`
PYTHON_ORG=`otool -L libsakurapy.so | grep Python | sed -e "s@^[ ${TAB}]*\(/.*/Python.framework/.*/Python\).*@\1@"`
PYTHON_NEW=`echo ${PYTHON_ORG} | sed -e "s@.*/\(Frameworks/Python.framework/.*\)@/Applications/CASA.app/Contents/\1@"`
echo ${PYTHON_ORG}
echo ${PYTHON_NEW}
install_name_tool -change ${PYTHON_ORG} ${PYTHON_NEW} libsakurapy.so

