#include <cstddef>
#include <map>
#include <vector>
#include <iostream>

#define PY_SSIZE_T_CLEAN
#include <Python.h>
// Author: Lajos Palanki -- Reused with permission from author
static std::map<std::vector<size_t>,bool> generated;

static std::vector<std::vector<size_t>> gen_bases(const size_t n, const size_t l, std::vector<size_t>* state = nullptr){
    std::vector<std::vector<size_t>> bases;
    std::vector<size_t> prev(l+1);
    if(!state){
        generated = std::map<std::vector<size_t>,bool>();
        if (l == 0)
            prev[0] = n;
        else {
                prev[0] = n-1;
                prev[l] = 1;
                
            }
    }else
        prev = *state;
    bases.push_back(prev);

    size_t j = 0;
    for(size_t i = l; i > 0; --i){
        j = 0;
        for(;j + 1 < i; ++j){
            if(prev[i] == 0 || prev[j] == 0) continue;
            std::vector<size_t> new_state = prev;
            --new_state[i],++new_state[i-1],--new_state[j],++new_state[j+1];
            if(generated[new_state]) continue;
            else generated[new_state] = true;
            for(const auto& elem : gen_bases(n, l, &new_state))
                bases.push_back(elem);
        }
    }

    return bases;
}

static PyObject* gen_bases_py_wrap(PyObject* self, PyObject* args){
    size_t l;
    size_t n;
    if(!PyArg_ParseTuple(args,"KK",&n,&l))
        return NULL;

    PyObject* py_bases = PyList_New(0);
    for (const auto& base : gen_bases(n,l)){
        PyObject* py_base = PyList_New(0);
        for (const auto& member : base)
            PyList_Append(py_base, PyLong_FromSize_t(member));
        PyList_Append(py_bases,py_base);
        Py_DECREF(py_base);
    }
    return py_bases;
}

static PyMethodDef constlmethods[] = {
    {"gen_bases",  gen_bases_py_wrap, METH_VARARGS,
     "Generate all bases having l angular moment with n particles"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef constlmodule = {
    PyModuleDef_HEAD_INIT,
    "const_l_gen",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    constlmethods
};

PyMODINIT_FUNC PyInit_const_l_gen(void){
    PyObject* m;

    m = PyModule_Create(&constlmodule);
    if(!m)
        return NULL;
    return m;
}

// int main(){
//     // for (size_t n = 1; n <= 10; ++n)
//     // for (size_t l = 0; l <= 90; ++l){
//         size_t n = 8;
//         size_t l = 76;
//         std::cout << n << " " << l << "\n";
//         std::vector<std::vector<size_t>> bases = gen_bases(n,l);
//         std::cout << bases.size() << "\n";
//     // }
//     return 0;
// }
