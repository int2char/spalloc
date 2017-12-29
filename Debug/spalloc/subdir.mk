################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../spalloc/main.cpp 

CU_SRCS += \
../spalloc/parallelor.cu 

CU_DEPS += \
./spalloc/parallelor.d 

OBJS += \
./spalloc/main.o \
./spalloc/parallelor.o 

CPP_DEPS += \
./spalloc/main.d 


# Each subdirectory must supply rules for building sources it contributes
spalloc/%.o: ../spalloc/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-8.0/bin/nvcc -I/home/zhang/allocate/spalloc/include -G -g -O0 -std=c++11 -gencode arch=compute_35,code=sm_35  -odir "spalloc" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-8.0/bin/nvcc -I/home/zhang/allocate/spalloc/include -G -g -O0 -std=c++11 --compile  -x c++ -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

spalloc/%.o: ../spalloc/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-8.0/bin/nvcc -I/home/zhang/allocate/spalloc/include -G -g -O0 -std=c++11 -gencode arch=compute_35,code=sm_35  -odir "spalloc" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-8.0/bin/nvcc -I/home/zhang/allocate/spalloc/include -G -g -O0 -std=c++11 --compile --relocatable-device-code=false -gencode arch=compute_35,code=compute_35 -gencode arch=compute_35,code=sm_35  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


