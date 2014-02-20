################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/cppad.cpp 

OBJS += \
./src/cppad.o 

CPP_DEPS += \
./src/cppad.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -DNDEBUG -I/home/kosari/catkin_ws/devel/include -I/opt/ros/hydro/include -I/home/kosari/Downloads/cppad-20140210 -I/home/kosari/Downloads/eigen -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


