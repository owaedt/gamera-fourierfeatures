
        
    
  #include "gameramodule.hpp"
  #include "knnmodule.hpp"

        #include "dft.hpp"
      #include "broken.hpp"
  
    #include <string>
  #include <stdexcept>
  #include "Python.h"
  #include <list>

  using namespace Gamera;
      using namespace FdToolkit;
  
        
    #ifndef _MSC_VER
  extern "C" {
    void init_broken(void);
  }
#endif

                static PyObject* call_fdbroken_a(PyObject* self, PyObject* args) {
            
      PyErr_Clear();
                                                                                                        Image* self_arg;
PyObject* self_pyarg;
      
                           int offset = -1;
         if (PyArg_ParseTuple(args,  "O|i:fdbroken_a",&self_pyarg, &offset) <= 0)
           return 0;
      
              if (!is_ImageObject(self_pyarg)) {
          PyErr_SetString(PyExc_TypeError, "Argument 'self' must be an image");
          return 0;
        }
        self_arg = ((Image*)((RectObject*)self_pyarg)->m_x);
        image_get_fv(self_pyarg, &self_arg->features, &self_arg->features_len);
              
               feature_t* feature_buffer = 0;
         if (offset < 0) {
           feature_buffer = new feature_t[ 60];
         } else {
           if (self_arg->features_len < offset + 60) {
             PyErr_Format(PyExc_ValueError, "Offset as given (%d) will cause data to be written outside of array of length (%d).  Perhaps the feature array is not initialised?", offset, (int)self_arg->features_len);
             return 0;
           }
           feature_buffer = self_arg->features + offset;
         }
         switch(get_image_combination(self_pyarg)) {
case ONEBITIMAGEVIEW:
fdbroken_a(*((OneBitImageView*)self_arg), feature_buffer);break;
case CC:
fdbroken_a(*((Cc*)self_arg), feature_buffer);break;
case ONEBITRLEIMAGEVIEW:
fdbroken_a(*((OneBitRleImageView*)self_arg), feature_buffer);break;
case RLECC:
fdbroken_a(*((RleCc*)self_arg), feature_buffer);break;
case MLCC:
fdbroken_a(*((MlCc*)self_arg), feature_buffer);break;
default:
PyErr_Format(PyExc_TypeError,"The 'self' argument of 'fdbroken_a' can not have pixel type '%s'. Acceptable values are ONEBIT, ONEBIT, ONEBIT, ONEBIT, and ONEBIT.", get_pixel_type_name(self_pyarg));
return 0;
}
      
               if (offset < 0) {
           PyObject* str = PyBytes_FromStringAndSize((char*)feature_buffer, 60* sizeof(feature_t));
           if (str != 0) {
                            PyObject* array_init = get_ArrayInit();
              if (array_init == 0)
                return 0;
              PyObject* array = PyObject_CallFunction(
                    array_init, (char *)"sO", (char *)"d", str);
              Py_XDECREF(str);
              delete[] feature_buffer;
              return array;
           } else {
             delete[] feature_buffer;
             return 0;
           }
         } else {
           Py_XINCREF(Py_None);
           return Py_None;
         }
            }
                static PyObject* call_fdbroken_a_phase(PyObject* self, PyObject* args) {
            
      PyErr_Clear();
                                                                                                        Image* self_arg;
PyObject* self_pyarg;
      
                           int offset = -1;
         if (PyArg_ParseTuple(args,  "O|i:fdbroken_a_phase",&self_pyarg, &offset) <= 0)
           return 0;
      
              if (!is_ImageObject(self_pyarg)) {
          PyErr_SetString(PyExc_TypeError, "Argument 'self' must be an image");
          return 0;
        }
        self_arg = ((Image*)((RectObject*)self_pyarg)->m_x);
        image_get_fv(self_pyarg, &self_arg->features, &self_arg->features_len);
              
               feature_t* feature_buffer = 0;
         if (offset < 0) {
           feature_buffer = new feature_t[ 60];
         } else {
           if (self_arg->features_len < offset + 60) {
             PyErr_Format(PyExc_ValueError, "Offset as given (%d) will cause data to be written outside of array of length (%d).  Perhaps the feature array is not initialised?", offset, (int)self_arg->features_len);
             return 0;
           }
           feature_buffer = self_arg->features + offset;
         }
         switch(get_image_combination(self_pyarg)) {
case ONEBITIMAGEVIEW:
fdbroken_a_phase(*((OneBitImageView*)self_arg), feature_buffer);break;
case CC:
fdbroken_a_phase(*((Cc*)self_arg), feature_buffer);break;
case ONEBITRLEIMAGEVIEW:
fdbroken_a_phase(*((OneBitRleImageView*)self_arg), feature_buffer);break;
case RLECC:
fdbroken_a_phase(*((RleCc*)self_arg), feature_buffer);break;
case MLCC:
fdbroken_a_phase(*((MlCc*)self_arg), feature_buffer);break;
default:
PyErr_Format(PyExc_TypeError,"The 'self' argument of 'fdbroken_a_phase' can not have pixel type '%s'. Acceptable values are ONEBIT, ONEBIT, ONEBIT, ONEBIT, and ONEBIT.", get_pixel_type_name(self_pyarg));
return 0;
}
      
               if (offset < 0) {
           PyObject* str = PyBytes_FromStringAndSize((char*)feature_buffer, 60* sizeof(feature_t));
           if (str != 0) {
                            PyObject* array_init = get_ArrayInit();
              if (array_init == 0)
                return 0;
              PyObject* array = PyObject_CallFunction(
                    array_init, (char *)"sO", (char *)"d", str);
              Py_XDECREF(str);
              delete[] feature_buffer;
              return array;
           } else {
             delete[] feature_buffer;
             return 0;
           }
         } else {
           Py_XINCREF(Py_None);
           return Py_None;
         }
            }
                static PyObject* call_fdbroken_a_phase_s1(PyObject* self, PyObject* args) {
            
      PyErr_Clear();
                                                                                                        Image* self_arg;
PyObject* self_pyarg;
      
                           int offset = -1;
         if (PyArg_ParseTuple(args,  "O|i:fdbroken_a_phase_s1",&self_pyarg, &offset) <= 0)
           return 0;
      
              if (!is_ImageObject(self_pyarg)) {
          PyErr_SetString(PyExc_TypeError, "Argument 'self' must be an image");
          return 0;
        }
        self_arg = ((Image*)((RectObject*)self_pyarg)->m_x);
        image_get_fv(self_pyarg, &self_arg->features, &self_arg->features_len);
              
               feature_t* feature_buffer = 0;
         if (offset < 0) {
           feature_buffer = new feature_t[ 60];
         } else {
           if (self_arg->features_len < offset + 60) {
             PyErr_Format(PyExc_ValueError, "Offset as given (%d) will cause data to be written outside of array of length (%d).  Perhaps the feature array is not initialised?", offset, (int)self_arg->features_len);
             return 0;
           }
           feature_buffer = self_arg->features + offset;
         }
         switch(get_image_combination(self_pyarg)) {
case ONEBITIMAGEVIEW:
fdbroken_a_phase_s1(*((OneBitImageView*)self_arg), feature_buffer);break;
case CC:
fdbroken_a_phase_s1(*((Cc*)self_arg), feature_buffer);break;
case ONEBITRLEIMAGEVIEW:
fdbroken_a_phase_s1(*((OneBitRleImageView*)self_arg), feature_buffer);break;
case RLECC:
fdbroken_a_phase_s1(*((RleCc*)self_arg), feature_buffer);break;
case MLCC:
fdbroken_a_phase_s1(*((MlCc*)self_arg), feature_buffer);break;
default:
PyErr_Format(PyExc_TypeError,"The 'self' argument of 'fdbroken_a_phase_s1' can not have pixel type '%s'. Acceptable values are ONEBIT, ONEBIT, ONEBIT, ONEBIT, and ONEBIT.", get_pixel_type_name(self_pyarg));
return 0;
}
      
               if (offset < 0) {
           PyObject* str = PyBytes_FromStringAndSize((char*)feature_buffer, 60* sizeof(feature_t));
           if (str != 0) {
                            PyObject* array_init = get_ArrayInit();
              if (array_init == 0)
                return 0;
              PyObject* array = PyObject_CallFunction(
                    array_init, (char *)"sO", (char *)"d", str);
              Py_XDECREF(str);
              delete[] feature_buffer;
              return array;
           } else {
             delete[] feature_buffer;
             return 0;
           }
         } else {
           Py_XINCREF(Py_None);
           return Py_None;
         }
            }
                static PyObject* call_fdbroken_b(PyObject* self, PyObject* args) {
            
      PyErr_Clear();
                                                                                                        Image* self_arg;
PyObject* self_pyarg;
      
                           int offset = -1;
         if (PyArg_ParseTuple(args,  "O|i:fdbroken_b",&self_pyarg, &offset) <= 0)
           return 0;
      
              if (!is_ImageObject(self_pyarg)) {
          PyErr_SetString(PyExc_TypeError, "Argument 'self' must be an image");
          return 0;
        }
        self_arg = ((Image*)((RectObject*)self_pyarg)->m_x);
        image_get_fv(self_pyarg, &self_arg->features, &self_arg->features_len);
              
               feature_t* feature_buffer = 0;
         if (offset < 0) {
           feature_buffer = new feature_t[ 60];
         } else {
           if (self_arg->features_len < offset + 60) {
             PyErr_Format(PyExc_ValueError, "Offset as given (%d) will cause data to be written outside of array of length (%d).  Perhaps the feature array is not initialised?", offset, (int)self_arg->features_len);
             return 0;
           }
           feature_buffer = self_arg->features + offset;
         }
         switch(get_image_combination(self_pyarg)) {
case ONEBITIMAGEVIEW:
fdbroken_b(*((OneBitImageView*)self_arg), feature_buffer);break;
case CC:
fdbroken_b(*((Cc*)self_arg), feature_buffer);break;
case ONEBITRLEIMAGEVIEW:
fdbroken_b(*((OneBitRleImageView*)self_arg), feature_buffer);break;
case RLECC:
fdbroken_b(*((RleCc*)self_arg), feature_buffer);break;
case MLCC:
fdbroken_b(*((MlCc*)self_arg), feature_buffer);break;
default:
PyErr_Format(PyExc_TypeError,"The 'self' argument of 'fdbroken_b' can not have pixel type '%s'. Acceptable values are ONEBIT, ONEBIT, ONEBIT, ONEBIT, and ONEBIT.", get_pixel_type_name(self_pyarg));
return 0;
}
      
               if (offset < 0) {
           PyObject* str = PyBytes_FromStringAndSize((char*)feature_buffer, 60* sizeof(feature_t));
           if (str != 0) {
                            PyObject* array_init = get_ArrayInit();
              if (array_init == 0)
                return 0;
              PyObject* array = PyObject_CallFunction(
                    array_init, (char *)"sO", (char *)"d", str);
              Py_XDECREF(str);
              delete[] feature_buffer;
              return array;
           } else {
             delete[] feature_buffer;
             return 0;
           }
         } else {
           Py_XINCREF(Py_None);
           return Py_None;
         }
            }
                static PyObject* call_fdbroken_c(PyObject* self, PyObject* args) {
            
      PyErr_Clear();
                                                                                                        Image* self_arg;
PyObject* self_pyarg;
      
                           int offset = -1;
         if (PyArg_ParseTuple(args,  "O|i:fdbroken_c",&self_pyarg, &offset) <= 0)
           return 0;
      
              if (!is_ImageObject(self_pyarg)) {
          PyErr_SetString(PyExc_TypeError, "Argument 'self' must be an image");
          return 0;
        }
        self_arg = ((Image*)((RectObject*)self_pyarg)->m_x);
        image_get_fv(self_pyarg, &self_arg->features, &self_arg->features_len);
              
               feature_t* feature_buffer = 0;
         if (offset < 0) {
           feature_buffer = new feature_t[ 60];
         } else {
           if (self_arg->features_len < offset + 60) {
             PyErr_Format(PyExc_ValueError, "Offset as given (%d) will cause data to be written outside of array of length (%d).  Perhaps the feature array is not initialised?", offset, (int)self_arg->features_len);
             return 0;
           }
           feature_buffer = self_arg->features + offset;
         }
         switch(get_image_combination(self_pyarg)) {
case ONEBITIMAGEVIEW:
fdbroken_c(*((OneBitImageView*)self_arg), feature_buffer);break;
case CC:
fdbroken_c(*((Cc*)self_arg), feature_buffer);break;
case ONEBITRLEIMAGEVIEW:
fdbroken_c(*((OneBitRleImageView*)self_arg), feature_buffer);break;
case RLECC:
fdbroken_c(*((RleCc*)self_arg), feature_buffer);break;
case MLCC:
fdbroken_c(*((MlCc*)self_arg), feature_buffer);break;
default:
PyErr_Format(PyExc_TypeError,"The 'self' argument of 'fdbroken_c' can not have pixel type '%s'. Acceptable values are ONEBIT, ONEBIT, ONEBIT, ONEBIT, and ONEBIT.", get_pixel_type_name(self_pyarg));
return 0;
}
      
               if (offset < 0) {
           PyObject* str = PyBytes_FromStringAndSize((char*)feature_buffer, 60* sizeof(feature_t));
           if (str != 0) {
                            PyObject* array_init = get_ArrayInit();
              if (array_init == 0)
                return 0;
              PyObject* array = PyObject_CallFunction(
                    array_init, (char *)"sO", (char *)"d", str);
              Py_XDECREF(str);
              delete[] feature_buffer;
              return array;
           } else {
             delete[] feature_buffer;
             return 0;
           }
         } else {
           Py_XINCREF(Py_None);
           return Py_None;
         }
            }
                static PyObject* call_fdbroken_c_phase(PyObject* self, PyObject* args) {
            
      PyErr_Clear();
                                                                                                        Image* self_arg;
PyObject* self_pyarg;
      
                           int offset = -1;
         if (PyArg_ParseTuple(args,  "O|i:fdbroken_c_phase",&self_pyarg, &offset) <= 0)
           return 0;
      
              if (!is_ImageObject(self_pyarg)) {
          PyErr_SetString(PyExc_TypeError, "Argument 'self' must be an image");
          return 0;
        }
        self_arg = ((Image*)((RectObject*)self_pyarg)->m_x);
        image_get_fv(self_pyarg, &self_arg->features, &self_arg->features_len);
              
               feature_t* feature_buffer = 0;
         if (offset < 0) {
           feature_buffer = new feature_t[ 60];
         } else {
           if (self_arg->features_len < offset + 60) {
             PyErr_Format(PyExc_ValueError, "Offset as given (%d) will cause data to be written outside of array of length (%d).  Perhaps the feature array is not initialised?", offset, (int)self_arg->features_len);
             return 0;
           }
           feature_buffer = self_arg->features + offset;
         }
         switch(get_image_combination(self_pyarg)) {
case ONEBITIMAGEVIEW:
fdbroken_c_phase(*((OneBitImageView*)self_arg), feature_buffer);break;
case CC:
fdbroken_c_phase(*((Cc*)self_arg), feature_buffer);break;
case ONEBITRLEIMAGEVIEW:
fdbroken_c_phase(*((OneBitRleImageView*)self_arg), feature_buffer);break;
case RLECC:
fdbroken_c_phase(*((RleCc*)self_arg), feature_buffer);break;
case MLCC:
fdbroken_c_phase(*((MlCc*)self_arg), feature_buffer);break;
default:
PyErr_Format(PyExc_TypeError,"The 'self' argument of 'fdbroken_c_phase' can not have pixel type '%s'. Acceptable values are ONEBIT, ONEBIT, ONEBIT, ONEBIT, and ONEBIT.", get_pixel_type_name(self_pyarg));
return 0;
}
      
               if (offset < 0) {
           PyObject* str = PyBytes_FromStringAndSize((char*)feature_buffer, 60* sizeof(feature_t));
           if (str != 0) {
                            PyObject* array_init = get_ArrayInit();
              if (array_init == 0)
                return 0;
              PyObject* array = PyObject_CallFunction(
                    array_init, (char *)"sO", (char *)"d", str);
              Py_XDECREF(str);
              delete[] feature_buffer;
              return array;
           } else {
             delete[] feature_buffer;
             return 0;
           }
         } else {
           Py_XINCREF(Py_None);
           return Py_None;
         }
            }
                static PyObject* call_fdbroken_c_phase_s1(PyObject* self, PyObject* args) {
            
      PyErr_Clear();
                                                                                                        Image* self_arg;
PyObject* self_pyarg;
      
                           int offset = -1;
         if (PyArg_ParseTuple(args,  "O|i:fdbroken_c_phase_s1",&self_pyarg, &offset) <= 0)
           return 0;
      
              if (!is_ImageObject(self_pyarg)) {
          PyErr_SetString(PyExc_TypeError, "Argument 'self' must be an image");
          return 0;
        }
        self_arg = ((Image*)((RectObject*)self_pyarg)->m_x);
        image_get_fv(self_pyarg, &self_arg->features, &self_arg->features_len);
              
               feature_t* feature_buffer = 0;
         if (offset < 0) {
           feature_buffer = new feature_t[ 60];
         } else {
           if (self_arg->features_len < offset + 60) {
             PyErr_Format(PyExc_ValueError, "Offset as given (%d) will cause data to be written outside of array of length (%d).  Perhaps the feature array is not initialised?", offset, (int)self_arg->features_len);
             return 0;
           }
           feature_buffer = self_arg->features + offset;
         }
         switch(get_image_combination(self_pyarg)) {
case ONEBITIMAGEVIEW:
fdbroken_c_phase_s1(*((OneBitImageView*)self_arg), feature_buffer);break;
case CC:
fdbroken_c_phase_s1(*((Cc*)self_arg), feature_buffer);break;
case ONEBITRLEIMAGEVIEW:
fdbroken_c_phase_s1(*((OneBitRleImageView*)self_arg), feature_buffer);break;
case RLECC:
fdbroken_c_phase_s1(*((RleCc*)self_arg), feature_buffer);break;
case MLCC:
fdbroken_c_phase_s1(*((MlCc*)self_arg), feature_buffer);break;
default:
PyErr_Format(PyExc_TypeError,"The 'self' argument of 'fdbroken_c_phase_s1' can not have pixel type '%s'. Acceptable values are ONEBIT, ONEBIT, ONEBIT, ONEBIT, and ONEBIT.", get_pixel_type_name(self_pyarg));
return 0;
}
      
               if (offset < 0) {
           PyObject* str = PyBytes_FromStringAndSize((char*)feature_buffer, 60* sizeof(feature_t));
           if (str != 0) {
                            PyObject* array_init = get_ArrayInit();
              if (array_init == 0)
                return 0;
              PyObject* array = PyObject_CallFunction(
                    array_init, (char *)"sO", (char *)"d", str);
              Py_XDECREF(str);
              delete[] feature_buffer;
              return array;
           } else {
             delete[] feature_buffer;
             return 0;
           }
         } else {
           Py_XINCREF(Py_None);
           return Py_None;
         }
            }
      
          static PyMethodDef _broken_methods[] = {
                  {  "fdbroken_a",
          call_fdbroken_a, METH_VARARGS,
           "**fdbroken_a** ()\n\nFourier descriptor for broken shapes according to Eq.(19) in\n[Dalitz2013]_. The coefficient c_0 is used for scale normalisation.\n\nThe absolute values of \nthe coefficients A(0), A(N-1) A(1), A(N-2), ... are returned.\n\nReference:\n\n  .. [Dalitz2013] Dalitz, C., Brandt, C., Goebbels, S., Kolanus, D.:\n     \"Fourier Descriptors for Broken Shapes.\" EURASIP Journal on\n     Advances in Signal Processing 2013:161, 2013\n\n.. note:: In this study, *fdbroken_a* was the best performing Fourier\n          descriptor on the NEUMES data set of broken shapes.\n"        },
                        {  "fdbroken_a_phase",
          call_fdbroken_a_phase, METH_VARARGS,
           "**fdbroken_a_phase** ()\n\nFourier descriptor for broken shapes according to Eq. (19) in\n[Dalitz2013]_. The coefficients c_r and c_s are used for scale and phase\nnormalisation. c_r and c_s are the two coefficients with the largest and\nsecond largest value (r >= 0, 0 < s < N/2).\n\nThe complex values of the coefficients A(0), A(N/2-1) A(1), A(N/2-2), ...\nare returned as alternating real and imaginary parts."        },
                        {  "fdbroken_a_phase_s1",
          call_fdbroken_a_phase_s1, METH_VARARGS,
           "**fdbroken_a_phase_s1** ()\n\nFourier descriptor for broken shapes according to Eq. (19) in\n[Dalitz2013]_. \nThe coefficients c_0 and c_1 are used for scale and phase normalisation. \n\nThe complex values of the coefficients A(0), A(N/2-1) A(1), A(N/2-2), ... \nare returned as alternating real and imaginary parts."        },
                        {  "fdbroken_b",
          call_fdbroken_b, METH_VARARGS,
           "**fdbroken_b** ()\n\nFourier descriptor for broken shapes according to Eq. (20) in\n[Dalitz2013]_.\n\nThe values of the coefficients B(1), B(2), ... B(N) are returned."        },
                        {  "fdbroken_c",
          call_fdbroken_c, METH_VARARGS,
           "**fdbroken_c** ()\n\nFourier descriptor for broken shapes according to Eq. (22) in\n[Dalitz2013]_. The coefficient c_0 is used for scale normalisation.\n\nThe absolute values of the coefficients C(0), ... C(N/2-1) are returned."        },
                        {  "fdbroken_c_phase",
          call_fdbroken_c_phase, METH_VARARGS,
           "**fdbroken_c_phase** ()\n\nFourier descriptor for broken shapes according to Eq. (22) in\n[Dalitz2013]_. The coefficients c_r and c_s are used for scale and phase\nnormalisation. c_r and c_s are the two coefficients with the largest and\nsecond largest absolute value (r >= 0, 0 < s < N/2).\n\nThe complex values of the coefficients C(0), \n... C(N/4-1) are returned as alternating real and imaginary parts."        },
                        {  "fdbroken_c_phase_s1",
          call_fdbroken_c_phase_s1, METH_VARARGS,
           "**fdbroken_c_phase_s1** ()\n\nFourier descriptor for broken shapes according to Eq. (22) in\n[Dalitz2013]_. The coefficients c_0 = 1/N  sum(r(t)) is used for scale\nnormalisation. The coefficient c_1 is used for phase normalization.\n\nThe complex values of the coefficients C(0), ... C(N/4-1) are returned as\nalternating real and imaginary parts."        },
              { nullptr }
  };
  
  static struct PyModuleDef module_brokenDef = {
        PyModuleDef_HEAD_INIT,
        "_broken",
        nullptr,
        -1,
        _broken_methods,
        nullptr,
        nullptr,
        nullptr,
        nullptr
  };


  PyMODINIT_FUNC PyInit__broken(void) {
    return PyModule_Create(&module_brokenDef);
  }
  

