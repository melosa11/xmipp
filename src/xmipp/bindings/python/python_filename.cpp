/***************************************************************************
 *
 * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "python_filename.h"
#include "core/xmipp_filename.h"
#include "core/xmipp_image_macros.h"
#include "core/xmipp_image_extension.h"
#include "core/xmipp_error.h"

/***************************************************************/
/*                            FileName                         */
/**************************************************************/

/* FileName methods */
PyMethodDef FileName_methods[] =
{
    { "compose", (PyCFunction) FileName_compose, METH_VARARGS,
      "Compose from root, number and extension OR prefix with number @" },
    { "composeBlock", (PyCFunction) FileName_composeBlock, METH_VARARGS,
      "Compose from blockname, number, root and extension" },
    { "decompose", (PyCFunction) FileName_decompose, METH_NOARGS,
      "Decompose filenames with @. Mainly from selfiles" },
    { "getBaseName", (PyCFunction) FileName_getBaseName, METH_NOARGS,
      "Get the base name from a FileName" },
    { "getExtension", (PyCFunction) FileName_getExtension, METH_NOARGS,
      "Get the last extension from a FileName" },
    { "getNumber", (PyCFunction) FileName_getNumber, METH_NOARGS,
      "Get the number from a FileName" },
    { "isInStack", (PyCFunction) FileName_isInStack, METH_NOARGS,
      "True if filename has stack format" },
    { "exists", (PyCFunction) FileName_exists, METH_NOARGS,
      "True if FileName exists" },
    { "isMetaData", (PyCFunction) FileName_isMetaData, METH_NOARGS,
      "True if is a MetaData" },
    { "isImage", (PyCFunction) FileName_isImage, METH_NOARGS,
      "True if is an image" },
    { "isStar1", (PyCFunction) FileName_isStar1, METH_NOARGS,
      "True if is a Star1" },
    { "withoutExtension", (PyCFunction) FileName_withoutExtension, METH_NOARGS,
      "return filename without extension" },
    { "removeBlockName", (PyCFunction) FileName_removeBlockName, METH_NOARGS,
      "return filename without block" },
    { nullptr } /* Sentinel */
};//FileName_methods

/*FileName Type */
PyTypeObject FileNameType =
{
   PyObject_HEAD_INIT(nullptr)
   "xmipp.FileName", /*tp_name*/
   sizeof(FileNameObject), /*tp_basicsize*/
   0, /*tp_itemsize*/
   (destructor)FileName_dealloc, /*tp_dealloc*/
   0, /*tp_print*/
   0, /*tp_getattr*/
   0, /*tp_setattr*/
   0, /*tp_compare*/
   FileName_repr, /*tp_repr*/
   0, /*tp_as_number*/
   0, /*tp_as_sequence*/
   0, /*tp_as_mapping*/
   0, /*tp_hash */
   0, /*tp_call*/
   0, /*tp_str*/
   0, /*tp_getattro*/
   0, /*tp_setattro*/
   0, /*tp_as_buffer*/
   Py_TPFLAGS_DEFAULT, /*tp_flags*/
   "Python wrapper to Xmipp FileName class",/* tp_doc */
   0, /* tp_traverse */
   0, /* tp_clear */
   0, /* tp_richcompare */
   0, /* tp_weaklistoffset */
   0, /* tp_iter */
   0, /* tp_iternext */
   FileName_methods, /* tp_methods */
   0, /* tp_members */
   0, /* tp_getset */
   0, /* tp_base */
   0, /* tp_dict */
   0, /* tp_descr_get */
   0, /* tp_descr_set */
   0, /* tp_dictoffset */
   0, /* tp_init */
   0, /* tp_alloc */
   FileName_new, /* tp_new */
};//FileNameType

/* Destructor */
void FileName_dealloc(FileNameObject* self)
{
    delete self->filename;
    Py_TYPE(self)->tp_free((PyObject*)self);
}

/* Constructor */
PyObject *
FileName_new(PyTypeObject *type, PyObject *args, PyObject *kwargs)
{
    FileNameObject *self;

    self = (FileNameObject*)type->tp_alloc(type, 0);
    if (self != nullptr)
    {

        PyObject *input = nullptr; 
        PyObject *pyStr = nullptr;
        PyObject *pyStr2 = nullptr;
        char str[1024] = "";
        char ext[1024] = "";
        int number = ALL_IMAGES;
        if (PyArg_ParseTuple(args, "|Ois", &input, &number, &ext))
            //|| PyArg_ParseTuple(args, "|Os", &input, &ext)) FIXME
        {
            pyStr = PyObject_Str(input);
            if (pyStr != nullptr)
              {
                   const char *strExcType =  PyUnicode_AsUTF8(pyStr);
                   strcpy(str, strExcType);
              }
        }
        if (number != ALL_IMAGES)
            self->filename = new FileName(str, number, ext);
        else
            self->filename = new FileName(str);

    }
    return (PyObject *)self;
}

/* String representation */
PyObject *
FileName_repr(PyObject * obj)
{
    auto *self = (FileNameObject*) obj;
    return PyUnicode_FromString(self->filename->c_str());
}

/* compose */
PyObject *
FileName_compose(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    auto *self = (FileNameObject*) obj;
    PyObject* pyStr2 = nullptr;

    if (self != nullptr)
    {
        PyObject *input = nullptr;
        PyObject *pyStr = nullptr;
        PyObject *input2 = nullptr;
        char str[1024] = "";
        char * ext = nullptr;
        char str2[1024] = "";
        int number = -1;
        size_t n = PyTuple_Size(args);
        //"kk000001.xmp"
        if (n == 3 && PyArg_ParseTuple(args, "Ois", &input, &number, &ext))
        {
            pyStr = PyObject_Str(input);
            if (pyStr != nullptr)
                strcpy(str, PyUnicode_AsUTF8(pyStr));
            self->filename->compose(str, number, ext);
        }
        else if (n == 2  && PyArg_ParseTuple(args, "OO", &input, &input2))
        {
            if( PyUnicode_Check( input ) )
            {
                //"jj@kk.xmp"
                pyStr  = PyObject_Str(input);
                pyStr2 = PyObject_Str(input2);
                if (pyStr != nullptr)
                	strcpy(str, PyUnicode_AsUTF8(pyStr));
                if (pyStr2 != nullptr)
                	strcpy(str2, PyUnicode_AsUTF8(pyStr2));
                self->filename->compose(str, str2);
            }
            else if ( PyLong_Check( input ) )
            {
                //"1@kk.xmp"
                number=PyLong_AsLong(input);
                pyStr2 = PyObject_Str(input2);
                strcpy(str2, PyUnicode_AsUTF8(pyStr2));
                self->filename->compose(number, str2);
            }
            else
                return nullptr;
        }
        Py_RETURN_NONE;//Return None(similar to void in C)
    }
    Py_RETURN_NONE;//Return None(similar to void in C)
}
/* composeBlock */
PyObject *
FileName_composeBlock(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    auto *self = (FileNameObject*) obj;

    if (self != nullptr)
    {
        char root[1024] = "";
        char ext[32] = "";
        char block[1024] ="";
        int number = 1;
        PyArg_ParseTuple(args, "sis|s", &block, &number, &root, &ext);
        self->filename->composeBlock(block, number, root, ext);
    }
    Py_RETURN_NONE;//Return None(similar to void in C)
}

/* exists */
PyObject *
FileName_exists(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    auto *self = (FileNameObject*) obj;

    if (self->filename->existsTrim())
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

/* isInStack */
PyObject *
FileName_isInStack(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    auto *self = (FileNameObject*) obj;

    if (self->filename->isInStack())
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

/* isMetadata */
PyObject *
FileName_isMetaData(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    auto *self = (FileNameObject*) obj;
    try
    {
        if(self->filename->isMetaData(false))
        {
            Py_RETURN_TRUE;
        }
        else
        {
            Py_RETURN_FALSE;
        }
    }
    catch (XmippError &xe)
    {
        PyErr_SetString(PyXmippError, xe.msg.c_str());
    }
    return nullptr;
}

/* isImage */
PyObject *
FileName_isImage(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    if (isImage(FileName_Value(obj)))
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

/* isStar1 */
PyObject *
FileName_isStar1(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    auto *self = (FileNameObject*) obj;

    if (self->filename->isStar1(false))
        Py_RETURN_TRUE;
    else
        Py_RETURN_FALSE;
}

PyObject *
FileName_getExtension(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    auto *self = (FileNameObject*) obj;

    return PyUnicode_FromString(self->filename->getExtension().c_str());
}

PyObject *
FileName_getNumber(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    auto *self = (FileNameObject*) obj;

    return PyLong_FromLong(self->filename->getNumber());
}

PyObject *
FileName_getBaseName(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    auto *self = (FileNameObject*) obj;

    return PyUnicode_FromString(self->filename->getBaseName().c_str());
}

PyObject *
FileName_decompose(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    auto *self = (FileNameObject*) obj;
    size_t no;
    String str;
    self->filename->decompose(no, str);
    return Py_BuildValue("is", no, str.c_str());
}

PyObject *
FileName_withoutExtension(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    auto *self = (FileNameObject*) obj;
    return PyUnicode_FromString(self->filename->withoutExtension().c_str());
}

PyObject *
FileName_removeBlockName(PyObject *obj, PyObject *args, PyObject *kwargs)
{
    auto *self = (FileNameObject*) obj;
    return PyUnicode_FromString(self->filename->removeBlockName().c_str());
}
