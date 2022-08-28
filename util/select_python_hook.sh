#find out if we can use system's python or need included

if command -v X_PYBIN_X > /dev/null ; then
retval=$(X_PYBIN_X -c "
import sysconfig as s;
from os import path as p;
v = s.get_config_vars();
print(p.join(v['LIBDIR'], v['LDLIBRARY']))")
fi
if [ -f "${retval}" ]; then
    echo "Using system python"
    export LD_PRELOAD=${retval}
else
    echo "Using bundled python"
    export PYTHONHOME=$(dirname "$0")X_PYROOT_X
    export LD_LIBRARY_PATH="$PWD/optional/python":${LD_LIBRARY_PATH}
fi
