AML_MiniProject_1

Exercise 2.A.1
        Gyroscope measures angular velocity
        Accelerometer measures linear acceleration
    --> what definition for the pitch angular velocity? S.42, Lecture1
        --> Can assume the amplitude is higher around this axis, z axis
        Yo axis, roll axis (foot towards inside or outside, laterally) & pitch axis (axis toward the side)
        When walking: the pitch axis exhibits the highest variability in amplitude
    
Exercise 2.A.2
        IC = initial contact 
        TC = terminal contact
        mid-stance events during all foot-flat
    --> How to determine the beginning of the FF? (FF_idx)

Exercise 2.A.3
    --> How to deterimne mid-stance? FF_idx is doing that?
        Mid-stance = mid of the FF --> mean between IC & TC: IC + (TC - IC)/2
    --> Problem of integer indices in
        "plot(data.imu.right.time(ICR_idx),  rightGyro_AlignedwithFoot(ICR_idx,3), 'rx');"

Exercise 2.A.4
    --> very small values: very rapid walk?
    
Exercise 2.B.1
    --> dimension problem, if i = 1
            - time_1gaitCycle(1:i+1): length of 2
            - leftGyro_AlignedwithFoot(FFL_idx(8):(FFL_idx(9)+i), 3): size of 510*1
 
 
 
Exercise 3.A.1
    --> why should we do the average value?
    
Exercise 3.A.2
    --> Anatomical frame: is [0, -1, 0] right? NO
    Accelerometer: measures in free fall, static: measures reaction from the ground
    => not -1 BUT +1 for g!!!!!!
    
Exercise 3.A.3
    --> How is computed the rotation matrix? I saw the function at the end but how was it thought or what does it mean?
    --> How to choose which vector is A and which one is B?
    
Exercise 3.A.4 
    --> How to apply the rotation matrix? Why should we transpose?
    --> What do we expect as a result? Put the y signal on the zero?
    
Exercise 3.A.5
    --> Should we calculate 3 different matrix of rotation that we can apply to the different columns of the data.accelstat?
    
Exercise 3.B.1
    --> Way of thinking: is it correct?
    
Exercise 3.B.2
    --> How is computed the rotation matrix? 1st column = x, 2nd = y, 3rd = z?
    --> Should I take the inverse? Because it means: from TF to GF?
    
Exercise 3.B.3
    --> Procedure to find AF_x,y,z, is it correct?
    
Exercise 3.C.1 & 2
    loop on each time frame BUT
    --> keep the 3 dimensions?
    --> how to normalize?
    
Exercise 3.C.3
    1. Define anatomical axis of the foot
    2. Estimate the angle alpha <--- formula
    --> How to define the anatomical axis of the foot? (which corresponds to u)
    --> How to define n? 
    
