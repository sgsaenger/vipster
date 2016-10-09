%MappedType Vipster::Vec
{
%TypeHeaderCode
#include <vec.h>
%End

%ConvertFromTypeCode
return Py_BuildValue("(fff)", sipCpp[0], sipCpp[1], sipCpp[2]);
%End

%ConvertToTypeCode
if(sipIsErr == NULL) return PyTuple_Check(sipPy);
*sipCppPtr = new std::array<float,3>;
PyArg_ParseTuple(sipPy, "fff", sipCppPtr[0], sipCppPtr[1], sipCppPtr[2]);
return sipGetState(sipTransferObj);
%End
};
