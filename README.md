# sparse

At the base on transpilated CSparse from C to Go by [c4go](https://github.com/Konstantin8105/c4go).



Result of:
```
└─▪ cat sparse.go | grep -v "func()" | grep -v "function" | grep "func" | wc -l
75

└─▪ cat sparse.go | grep -v "func()" | grep -v "function" | grep "func" |  grep -v "//" | wc -l
13
```

#### **CSparse base on SuiteSparse 5.3.0**

http://faculty.cse.tamu.edu/davis/SuiteSparse/
