# Monitor GPU activity
``` 
watch -n 1 nvidia-smi
```

# Compile OpenACC code

## Matrix Multiplication
``` 
pgcc -o mat_multiplyGPU mat_multiplyGPU.c -acc -Minfo=acc -fast 
``` 

## Vector
``` 
pgcc vecadd.c -o vecaddACC -acc -Minfo=acc  -fast 
``` 

## ignoring openacc directives:
```
gcc vecadd.c -o vecadd -lm
```
