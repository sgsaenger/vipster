%Module vipster

%Include vec.sip

%MappedType std::string
{
%TypeHeaderCode
#include <string>
%End
 
%ConvertFromTypeCode
    // convert a std::string to a Python (unicode) string
    return PyUnicode_FromString(sipCpp->c_str());
%End
 
%ConvertToTypeCode
    // Allow a Python string whenever a string is
    // expected.
    if (sipIsErr == NULL) return PyUnicode_Check(sipPy);
    if (sipPy == Py_None) {
        *sipCppPtr = new std::string;
        return 1;
    }
    if (PyUnicode_Check(sipPy)) {
        *sipCppPtr = new std::string(PyUnicode_AsUTF8(sipPy));
        return 1;
    }
    return 0;
%End
};
