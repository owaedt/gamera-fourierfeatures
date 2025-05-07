
        
    
  #include "gameramodule.hpp"
  #include "knnmodule.hpp"

        #include "dft.hpp"
      #include "single.hpp"
  
    #include <string>
  #include <stdexcept>
  #include "Python.h"
  #include <list>

  using namespace Gamera;
      using namespace FdToolkit;
  
        
    #ifndef _MSC_VER
  extern "C" {
    void init_single(void);
  }
#endif

                static PyObject* call_fdsingle_complex_position(PyObject* self, PyObject* args) {
            
      PyErr_Clear();
                                                                                                        Image* self_arg;
PyObject* self_pyarg;
      
                           int offset = -1;
         if (PyArg_ParseTuple(args,  "O|i:fdsingle_complex_position",&self_pyarg, &offset) <= 0)
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
fdsingle_complex_position(*((OneBitImageView*)self_arg), feature_buffer);break;
case CC:
fdsingle_complex_position(*((Cc*)self_arg), feature_buffer);break;
case ONEBITRLEIMAGEVIEW:
fdsingle_complex_position(*((OneBitRleImageView*)self_arg), feature_buffer);break;
case RLECC:
fdsingle_complex_position(*((RleCc*)self_arg), feature_buffer);break;
case MLCC:
fdsingle_complex_position(*((MlCc*)self_arg), feature_buffer);break;
default:
PyErr_Format(PyExc_TypeError,"The 'self' argument of 'fdsingle_complex_position' can not have pixel type '%s'. Acceptable values are ONEBIT, ONEBIT, ONEBIT, ONEBIT, and ONEBIT.", get_pixel_type_name(self_pyarg));
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
                static PyObject* call_fdsingle_complex_position_r1(PyObject* self, PyObject* args) {
            
      PyErr_Clear();
                                                                                                        Image* self_arg;
PyObject* self_pyarg;
      
                           int offset = -1;
         if (PyArg_ParseTuple(args,  "O|i:fdsingle_complex_position_r1",&self_pyarg, &offset) <= 0)
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
fdsingle_complex_position_r1(*((OneBitImageView*)self_arg), feature_buffer);break;
case CC:
fdsingle_complex_position_r1(*((Cc*)self_arg), feature_buffer);break;
case ONEBITRLEIMAGEVIEW:
fdsingle_complex_position_r1(*((OneBitRleImageView*)self_arg), feature_buffer);break;
case RLECC:
fdsingle_complex_position_r1(*((RleCc*)self_arg), feature_buffer);break;
case MLCC:
fdsingle_complex_position_r1(*((MlCc*)self_arg), feature_buffer);break;
default:
PyErr_Format(PyExc_TypeError,"The 'self' argument of 'fdsingle_complex_position_r1' can not have pixel type '%s'. Acceptable values are ONEBIT, ONEBIT, ONEBIT, ONEBIT, and ONEBIT.", get_pixel_type_name(self_pyarg));
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
                static PyObject* call_fdsingle_complex_position_phase(PyObject* self, PyObject* args) {
            
      PyErr_Clear();
                                                                                                        Image* self_arg;
PyObject* self_pyarg;
      
                           int offset = -1;
         if (PyArg_ParseTuple(args,  "O|i:fdsingle_complex_position_phase",&self_pyarg, &offset) <= 0)
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
fdsingle_complex_position_phase(*((OneBitImageView*)self_arg), feature_buffer);break;
case CC:
fdsingle_complex_position_phase(*((Cc*)self_arg), feature_buffer);break;
case ONEBITRLEIMAGEVIEW:
fdsingle_complex_position_phase(*((OneBitRleImageView*)self_arg), feature_buffer);break;
case RLECC:
fdsingle_complex_position_phase(*((RleCc*)self_arg), feature_buffer);break;
case MLCC:
fdsingle_complex_position_phase(*((MlCc*)self_arg), feature_buffer);break;
default:
PyErr_Format(PyExc_TypeError,"The 'self' argument of 'fdsingle_complex_position_phase' can not have pixel type '%s'. Acceptable values are ONEBIT, ONEBIT, ONEBIT, ONEBIT, and ONEBIT.", get_pixel_type_name(self_pyarg));
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
                static PyObject* call_fdsingle_granlund(PyObject* self, PyObject* args) {
            
      PyErr_Clear();
                                                                                                        Image* self_arg;
PyObject* self_pyarg;
      
                           int offset = -1;
         if (PyArg_ParseTuple(args,  "O|i:fdsingle_granlund",&self_pyarg, &offset) <= 0)
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
fdsingle_granlund(*((OneBitImageView*)self_arg), feature_buffer);break;
case CC:
fdsingle_granlund(*((Cc*)self_arg), feature_buffer);break;
case ONEBITRLEIMAGEVIEW:
fdsingle_granlund(*((OneBitRleImageView*)self_arg), feature_buffer);break;
case RLECC:
fdsingle_granlund(*((RleCc*)self_arg), feature_buffer);break;
case MLCC:
fdsingle_granlund(*((MlCc*)self_arg), feature_buffer);break;
default:
PyErr_Format(PyExc_TypeError,"The 'self' argument of 'fdsingle_granlund' can not have pixel type '%s'. Acceptable values are ONEBIT, ONEBIT, ONEBIT, ONEBIT, and ONEBIT.", get_pixel_type_name(self_pyarg));
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
                static PyObject* call_fdsingle_complex_position_phase_r1(PyObject* self, PyObject* args) {
            
      PyErr_Clear();
                                                                                                        Image* self_arg;
PyObject* self_pyarg;
      
                           int offset = -1;
         if (PyArg_ParseTuple(args,  "O|i:fdsingle_complex_position_phase_r1",&self_pyarg, &offset) <= 0)
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
fdsingle_complex_position_phase_r1(*((OneBitImageView*)self_arg), feature_buffer);break;
case CC:
fdsingle_complex_position_phase_r1(*((Cc*)self_arg), feature_buffer);break;
case ONEBITRLEIMAGEVIEW:
fdsingle_complex_position_phase_r1(*((OneBitRleImageView*)self_arg), feature_buffer);break;
case RLECC:
fdsingle_complex_position_phase_r1(*((RleCc*)self_arg), feature_buffer);break;
case MLCC:
fdsingle_complex_position_phase_r1(*((MlCc*)self_arg), feature_buffer);break;
default:
PyErr_Format(PyExc_TypeError,"The 'self' argument of 'fdsingle_complex_position_phase_r1' can not have pixel type '%s'. Acceptable values are ONEBIT, ONEBIT, ONEBIT, ONEBIT, and ONEBIT.", get_pixel_type_name(self_pyarg));
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
                static PyObject* call_fdsingle_elliptic(PyObject* self, PyObject* args) {
            
      PyErr_Clear();
                                                                                                        Image* self_arg;
PyObject* self_pyarg;
      
                           int offset = -1;
         if (PyArg_ParseTuple(args,  "O|i:fdsingle_elliptic",&self_pyarg, &offset) <= 0)
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
fdsingle_elliptic(*((OneBitImageView*)self_arg), feature_buffer);break;
case CC:
fdsingle_elliptic(*((Cc*)self_arg), feature_buffer);break;
case ONEBITRLEIMAGEVIEW:
fdsingle_elliptic(*((OneBitRleImageView*)self_arg), feature_buffer);break;
case RLECC:
fdsingle_elliptic(*((RleCc*)self_arg), feature_buffer);break;
case MLCC:
fdsingle_elliptic(*((MlCc*)self_arg), feature_buffer);break;
default:
PyErr_Format(PyExc_TypeError,"The 'self' argument of 'fdsingle_elliptic' can not have pixel type '%s'. Acceptable values are ONEBIT, ONEBIT, ONEBIT, ONEBIT, and ONEBIT.", get_pixel_type_name(self_pyarg));
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
                static PyObject* call_fdsingle_centroid_distance(PyObject* self, PyObject* args) {
            
      PyErr_Clear();
                                                                                                        Image* self_arg;
PyObject* self_pyarg;
      
                           int offset = -1;
         if (PyArg_ParseTuple(args,  "O|i:fdsingle_centroid_distance",&self_pyarg, &offset) <= 0)
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
fdsingle_centroid_distance(*((OneBitImageView*)self_arg), feature_buffer);break;
case CC:
fdsingle_centroid_distance(*((Cc*)self_arg), feature_buffer);break;
case ONEBITRLEIMAGEVIEW:
fdsingle_centroid_distance(*((OneBitRleImageView*)self_arg), feature_buffer);break;
case RLECC:
fdsingle_centroid_distance(*((RleCc*)self_arg), feature_buffer);break;
case MLCC:
fdsingle_centroid_distance(*((MlCc*)self_arg), feature_buffer);break;
default:
PyErr_Format(PyExc_TypeError,"The 'self' argument of 'fdsingle_centroid_distance' can not have pixel type '%s'. Acceptable values are ONEBIT, ONEBIT, ONEBIT, ONEBIT, and ONEBIT.", get_pixel_type_name(self_pyarg));
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
                static PyObject* call_fdsingle_centroid_distance_phase(PyObject* self, PyObject* args) {
            
      PyErr_Clear();
                                                                                                        Image* self_arg;
PyObject* self_pyarg;
      
                           int offset = -1;
         if (PyArg_ParseTuple(args,  "O|i:fdsingle_centroid_distance_phase",&self_pyarg, &offset) <= 0)
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
fdsingle_centroid_distance_phase(*((OneBitImageView*)self_arg), feature_buffer);break;
case CC:
fdsingle_centroid_distance_phase(*((Cc*)self_arg), feature_buffer);break;
case ONEBITRLEIMAGEVIEW:
fdsingle_centroid_distance_phase(*((OneBitRleImageView*)self_arg), feature_buffer);break;
case RLECC:
fdsingle_centroid_distance_phase(*((RleCc*)self_arg), feature_buffer);break;
case MLCC:
fdsingle_centroid_distance_phase(*((MlCc*)self_arg), feature_buffer);break;
default:
PyErr_Format(PyExc_TypeError,"The 'self' argument of 'fdsingle_centroid_distance_phase' can not have pixel type '%s'. Acceptable values are ONEBIT, ONEBIT, ONEBIT, ONEBIT, and ONEBIT.", get_pixel_type_name(self_pyarg));
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
                static PyObject* call_fdsingle_real_position(PyObject* self, PyObject* args) {
            
      PyErr_Clear();
                                                                                                        Image* self_arg;
PyObject* self_pyarg;
      
                           int offset = -1;
         if (PyArg_ParseTuple(args,  "O|i:fdsingle_real_position",&self_pyarg, &offset) <= 0)
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
fdsingle_real_position(*((OneBitImageView*)self_arg), feature_buffer);break;
case CC:
fdsingle_real_position(*((Cc*)self_arg), feature_buffer);break;
case ONEBITRLEIMAGEVIEW:
fdsingle_real_position(*((OneBitRleImageView*)self_arg), feature_buffer);break;
case RLECC:
fdsingle_real_position(*((RleCc*)self_arg), feature_buffer);break;
case MLCC:
fdsingle_real_position(*((MlCc*)self_arg), feature_buffer);break;
default:
PyErr_Format(PyExc_TypeError,"The 'self' argument of 'fdsingle_real_position' can not have pixel type '%s'. Acceptable values are ONEBIT, ONEBIT, ONEBIT, ONEBIT, and ONEBIT.", get_pixel_type_name(self_pyarg));
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
      
          static PyMethodDef _single_methods[] = {
                  {  "fdsingle_complex_position",
          call_fdsingle_complex_position, METH_VARARGS,
           "**fdsingle_complex_position** ()\n\nComplex position Fourier descriptor for unbroken shapes according to \nEq. (6) in [Dalitz2013]_. The coefficient c_r is used for scale normalisation.\nc_r is the coefficient with the largest absolute value.\n\nThe absolute values\nof the coefficients l(1), l(N-1), l(2), l(N-2), ... are returned.\n\n.. note:: In the study [Dalitz2013]_, this was the best performing Fourier\n          descriptor on the MPEG7 data set of connected shapes.\n"        },
                        {  "fdsingle_complex_position_r1",
          call_fdsingle_complex_position_r1, METH_VARARGS,
           "**fdsingle_complex_position_r1** ()\n\nComplex position Fourier descriptor for unbroken shapes according to \nEq. (6) in [Dalitz2013]_. The coefficient c_1 is used for scale normalisation.\n\nThe absolute values of the coefficients l(1), l(N-1), l(2), l(N-2), ... are \nreturned."        },
                        {  "fdsingle_complex_position_phase",
          call_fdsingle_complex_position_phase, METH_VARARGS,
           "**fdsingle_complex_position_phase** ()\n\nComplex position Fourier descriptor for unbroken shapes according to \nEq. (6) in [Dalitz2013]_. The coefficients c_r and c_s are used for scale and\nphase normalisation. c_r and c_s are the two coefficients with the largest and \nsecond largest absolute value. (r >= 0, 0 < s < N/2).\n\nThe complex values of the coefficients l(1), l(N-1) l(2), l(N-2), ... are\nreturned as alternating real and imaginary parts."        },
                        {  "fdsingle_granlund",
          call_fdsingle_granlund, METH_VARARGS,
           "**fdsingle_granlund** ()\n\nGranlund's Fourier descriptor for unbroken shapes according to \nEq. (7) in [Dalitz2013]_. The coefficient c_1 is used for normalisation.\nc_r is the coefficient with the largest absolute value.\n\nThe complex values of the coefficients d(2), d(3), ... are returned as\nalternating real and imaginary parts."        },
                        {  "fdsingle_complex_position_phase_r1",
          call_fdsingle_complex_position_phase_r1, METH_VARARGS,
           "**fdsingle_complex_position_phase_r1** ()\n\nComplex position Fourier descriptor for unbroken shapes according to \nEq. (6) in [Dalitz2013]_. The coefficients c_1 and c_r are used for scale and\nphase normalisation. c_r is the coefficient with the largest absolute value.\n\nThe complex values of the coefficients l(1), l(N-1) l(2), l(N-2), ... are \nreturned as alternating real and imaginary parts."        },
                        {  "fdsingle_elliptic",
          call_fdsingle_elliptic, METH_VARARGS,
           "**fdsingle_elliptic** ()\n\nElliptic Fourier descriptor for unbroken shapes according to Eq. \n(12) in [Dalitz2013]_. The coefficient I_1 is used for normalisation.\n\nThe values of the coefficients I(1), J(1), k1(1), I(2), J(2), K(2), ... are\nreturned."        },
                        {  "fdsingle_centroid_distance",
          call_fdsingle_centroid_distance, METH_VARARGS,
           "**fdsingle_centroid_distance** ()\n\nCentroid distance Fourier descriptor for unbroken shapes according\nto Eq. (16) in [Dalitz2013]_. The coefficient c_0 is used for scale\nnormalisation.\n\nThe absolute values of the coefficients R(1), R(2), ... are returned."        },
                        {  "fdsingle_centroid_distance_phase",
          call_fdsingle_centroid_distance_phase, METH_VARARGS,
           "**fdsingle_centroid_distance_phase** ()\n\nCentroid distance Fourier descriptor for unbroken shapes according\nto Eq. (16) in [Dalitz2013]_. The coefficients c_0 and c_s are used for scale\nand phase normalisation. c_s is the coefficient with the largest absolute\nvalue (s >= 1).\n\nThe complex values of the coefficients R(1), R(2), ... are returned \nas alternating real and imaginary parts."        },
                        {  "fdsingle_real_position",
          call_fdsingle_real_position, METH_VARARGS,
           "**fdsingle_real_position** ()\n\nReal position Fourier descriptor for unbroken shapes according to \nEq. (13) in [Dalitz2013]_. The coefficient c_1 is used for normalisation.\n\nThe values of the coefficients with indexes starting at 1 are returned."        },
              { nullptr }
  };
  
  static struct PyModuleDef module_singleDef = {
        PyModuleDef_HEAD_INIT,
        "_single",
        nullptr,
        -1,
        _single_methods,
        nullptr,
        nullptr,
        nullptr,
        nullptr
  };


  PyMODINIT_FUNC PyInit__single(void) {
    return PyModule_Create(&module_singleDef);
  }
  

