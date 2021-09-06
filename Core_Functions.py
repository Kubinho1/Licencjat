
import topoly as tp


def file_to_list(file):
    '''
    Converts .xyz file to a list for KnotCore() and PrecisionCore() functions
    '''
    f = open(file)

    g = f.readlines()
    f.close()
    l = []
    k = []

    for i in g:
        x = i.split()
        if len(x) == 1:
            k.append(l)
            l = []
        else:    
            l.append(x)

    k.append(l)
    if k[0] == []:
        k.pop(0)
    return k




def KnotCore(file, margin = 10, delete1 = 1, delete2 = 1, rollback = 3, atoms_to_check = 2, closure_type = 2, knot_cutoff = 0.42, alex_tries = 200, max_cross_param = 50):      
    '''
    Takes a list of coordinates of a protein chain and returns:
    0, if there's an unknot
    1, if the type of a knot is unknown
    list containing the knot type of a protein chain and it's range ex. ['3_1', 30, 20] means a trefoil knot with a knot core ranging from the 31st atom to the 20th atom from the end.
    '''




    #file - input list from file_to_list() function
    #margin - how many atoms we are cutting at a given time
    #delete1 - starting index for cutting (ex. if we don't want to start from the beginning we can put delete1 = 49 and that starts the procedure from the 50th atom)
    #delete2 - same as above except it's counted from the back (delete2 = 30 means we start from the 30th atom counting from the back)
    #rollback - when we hit the point where there are no more knots, we can go back a set amount diagonally 
    #atoms_to_check - set a number of atoms to check after finding an unknot, so that it confirms there are really no more knots after our current point
    #closure_type - sets the 'closure' parameter for the tp.alexander() function
    #knot_cutoff - set the threshold value of what we count as a knot, ex. the value returned by alexander needs to be higher than our cutoff to be recognized as a knot (only on closure_type = 2)
    #alex_tries - sets the 'tries' parameter for tp.alexander() function
    #max_cross - sets the 'max_cross' parameter for tp.alexander() function

    rollback = rollback + margin   #hits an unknot, goes back to the last knot, then rolls back a set amount
    atoms = file                                     
    for i in range(len(atoms)):    #converts items of list to be floats           
            for j in range(3):
                atoms[i][j] = float(atoms[i][j])
    

    if closure_type == 2:             #case when chain is closed probabilistically

        knot_dict = tp.alexander(atoms, closure = closure_type, max_cross = max_cross_param, tries = alex_tries) #create a dict of knots and their probabilities
        try:                                                                                                     #check for the knot with the biggest value
            greatest_knot_value = max(knot_dict.values())
        except:
            greatest_knot_value = 0.0
  
        if greatest_knot_value >= knot_cutoff:                                                                   #check if knot probability is high enough

            for k in knot_dict:
            
                if knot_dict[k] == greatest_knot_value:
                  
                    knot_type = k                                                                                #set the knot type we'll be working with
        else:
        
            knot_type = '0_1'


        if knot_type == '0_1':                                                                                   #if unknot or unknown then end 

            return 0
         
        elif 'Unknown' in knot_type or 'TooManyCrossings' in knot_type:
        
            return 1

        else:                                                                                                    #We'll be cutting the chain from the left, checking the type of knot, 
                                                                                                                 #then cutting from the right, checking the type of knot and so on
            chain_len = len(atoms)                                                                             
            left_side_type = knot_type                                             
            right_side_type = knot_type                                             
            now_cutting_from_left = True 

        while left_side_type == knot_type and right_side_type == knot_type:                                      
  
            if now_cutting_from_left:                                                                            #left side cut

                if delete1 + delete2 > chain_len - 3:                                                            #at no point can the indices get so close that the resulting chain has length less than 3

                    return [knot_type, delete1, delete2] 

                cut = atoms[delete1:-delete2]                                                                    #cut the chain

                left_side_type_dict = tp.alexander(cut, closure = closure_type, max_cross = max_cross_param, tries = alex_tries)     #calculate knot probability
                
                try:                                                                                            
                    value_of_knot = float(left_side_type_dict[knot_type])
                except:
                    value_of_knot = 0.0

                if value_of_knot >= knot_cutoff:                                                                #if knot made the cutoff then continue from the right side
               
                    delete2 += margin
                    now_cutting_from_left = False

                else:                                                                                           #else check for knots deeper in the chain (how many times to check - atoms_to_check)
                    left_side_type = None
                    at = 1
 
                    while at <= atoms_to_check and value_of_knot < knot_cutoff:

                        cut = atoms[delete1+at*margin:-delete2-at*margin]

                        if len(cut) == 0:

                            break

                        else:

                            try:
                                value_of_knot = float(tp.alexander(cut, closure = closure_type, max_cross = max_cross_param, tries = alex_tries)[knot_type])
                            except:
                                value_of_knot = 0.0
 
                            at += 1

                    if value_of_knot >= knot_cutoff:                                                          #if there are more knots down the line, then continue cutting from the right
 
                        left_side_type = knot_type
                        delete2 += margin
                        now_cutting_from_left = False

                    else:                                                                                     #else roll back and go out of the loop 

                        delete1 -= rollback
                        delete2 -= rollback

                        if delete1 < 1:                                                                       #make sure the indices aren't negative
                            delete1 = 1
                        if delete2 < 1:
                            delete2 = 1
                     

            else:                                                                                             #The same exact procedure as above except we cut from the end of the chain
              
                if delete1 + delete2 > chain_len - 3:

                    return [knot_type, delete1, delete2] 

                cut = atoms[delete1:-delete2]

                right_side_type_dict = tp.alexander(cut, closure = closure_type, max_cross = max_cross_param, tries = alex_tries)
                
                try:
                    value_of_knot = float(right_side_type_dict[knot_type])
                except:
                    value_of_knot = 0.0

                if value_of_knot >= knot_cutoff:
               
                    delete1 += margin
                    now_cutting_from_left = True

                else:
                    right_side_type = None
                    at = 1
                    
                    while at <= atoms_to_check and value_of_knot < knot_cutoff:

                        cut = atoms[delete1+at*margin:-delete2-at*margin]

                        if len(cut) == 0:

                            break

                        else:

                            try:
                                value_of_knot = float(tp.alexander(cut, closure = closure_type, max_cross = max_cross_param, tries = alex_tries)[knot_type])
                            except:
                                value_of_knot = 0.0
                            
                            at += 1

                    if value_of_knot >= knot_cutoff:
                        
                        right_side_type = knot_type
                        delete1 += margin
                        now_cutting_from_left = True

                    else:

                        delete1 -= rollback
                        delete2 -= rollback

                        if delete1 < 1:
                            delete1 = 1
                        if delete2 < 1:
                            delete2 = 1
                     

        left_side_type = knot_type                  
        right_side_type = knot_type
        temp_val = delete1                         #in the following loop the 'delete1' variable gets changed so we need to store it for later

        while left_side_type == knot_type:         #The same procedure as before, except we cut exclusively from the left. The resulting value of the left index will be our final value
            
            if delete1 + delete2 > chain_len - 3:
          
                return [typ, delete1, delete2]

            cut = atoms[delete1:-delete2]

            left_side_type_dict = tp.alexander(cut, closure = closure_type, max_cross = max_cross_param, tries = alex_tries)
           
            try:
                value_of_knot = float(left_side_type_dict[knot_type])
            except:
                value_of_knot = 0.0

            if value_of_knot >= knot_cutoff:

                delete1 += margin

            else:
                left_side_type = None
                at = 1
                
                while at <= atoms_to_check and value_of_knot < knot_cutoff:

                    cut = atoms[delete1+at*margin:-delete2]
 
                    if len(cut) == 0:

                        break

                    else:

                        try:
                            value_of_knot = float(tp.alexander(cut, closure = closure_type, max_cross = max_cross_param, tries = alex_tries)[knot_type])
                        except:
                            value_of_knot = 0.0
                       
                        at += 1
                        
                    if value_of_knot >= knot_cutoff:
                       
                        left_side_type = knot_type
                        delete1 += margin

                    
               
                delete1 -= margin
                final_left = delete1             #This is the left border of the knot's core

        delete1 = temp_val

        while right_side_type == knot_type:               #The same procedure as before, except we cut exclusively from the right. The resulting value of the right index will be our final value
            
            if delete1 + delete2 > chain_len - 3:
          
                return [typ, delete1, delete2]

            cut = atoms[delete1:-delete2]

            right_side_type_dict = tp.alexander(cut, closure = closure_type, max_cross = max_cross_param, tries = alex_tries)
          
            try:
                value_of_knot = float(right_side_type_dict[knot_type])
            except:
                value_of_knot = 0.0

            if value_of_knot >= knot_cutoff:

                delete2 += margin

            else:
                right_side_type = None
                at = 1
                
                while at <= atoms_to_check and value_of_knot < knot_cutoff:

                    cut = atoms[delete1:-delete2-at*margin]
 
                    if len(cut) == 0:

                        break

                    else:

                        try:
                            value_of_knot = float(tp.alexander(cut, closure = closure_type, max_cross = max_cross_param, tries = alex_tries)[knot_type])
                        except:
                            value_of_knot = 0.0
                        
                        at += 1
 
                    if value_of_knot >= knot_cutoff:
                        
                        right_side_type = knot_type
                        delete2 += margin

                    
                
                delete2 -= margin
                final_right = delete2              
        if final_left < 0:
            final_left = 0
        if final_right < 1:
            final_right = 1
        return [knot_type, final_left, final_right]

    elif closure_type == 1 or closure_type == 0:                                              #Case where the chain is closed via the deterministic method, everything works the same except we dont't make 
                                                                                              #dictionaries of probabilities and instead use one the single return value of alexander
        knot_type = tp.alexander(atoms, closure = closure_type, max_cross = max_cross_param)
 
        if knot_type == '0_1':

            return 0
         
        elif 'Unknown' in knot_type or 'TooManyCrossings' in knot_type:
        
            return 1

        else:

            chain_len = len(atoms)
            left_side_type = knot_type                                             
            right_side_type = knot_type                                             
            now_cutting_from_left = True

        while left_side_type == knot_type and right_side_type == knot_type: 
            
            if now_cutting_from_left:

                if delete1 + delete2 > chain_len - 3:

                    return [knot_type, delete1, delete2] 

                cut = atoms[delete1:-delete2]

                left_side_type = tp.alexander(cut, closure = closure_type, max_cross = max_cross_param)
                

                if left_side_type == knot_type:
               
                    delete2 += margin
                    now_cutting_from_left = False

                else:
                    
                    at = 1
                    
                    while at <= atoms_to_check and left_side_type != knot_type:

                        cut = atoms[delete1+at*margin:-delete2-at*margin]
                        
                        if len(cut) == 0:

                            break

                        else:
                            
                            left_side_type = tp.alexander(cut, closure = closure_type, max_cross = max_cross_param)
                           
                            at += 1

                    if left_side_type == knot_type:
                      
                        delete2 += margin
                        now_cutting_from_left = False

                    else:

                        delete1 -= rollback
                        delete2 -= rollback

                        if delete1 < 1:
                            delete1 = 1
                        if delete2 < 1:
                            delete2 = 1
                     

            else:
                
                if delete1 + delete2 > chain_len - 3:

                    return [knot_type, delete1, delete2] 

                cut = atoms[delete1:-delete2]
                
                right_side_type = tp.alexander(cut, closure = closure_type, max_cross = max_cross_param)
                
                if right_side_type == knot_type:
               
                    delete1 += margin
                    now_cutting_from_left = True


                else:
                    
                    at = 1
                    
                    while at <= atoms_to_check and right_side_type != knot_type:

                        cut = atoms[delete1+at*margin:-delete2-at*margin]
                        
                        if len(cut) == 0:

                            break

                        else:

                            right_side_type = tp.alexander(cut, closure = closure_type, max_cross = max_cross_param)
                          
                            at += 1

                    if right_side_type == knot_type:
                       
                        delete1 += margin
                        now_cutting_from_left = True

                    else:

                        delete1 -= rollback
                        delete2 -= rollback

                        if delete1 < 1:
                            delete1 = 1
                        if delete2 < 1:
                            delete2 = 1
        
        left_side_type = knot_type
        right_side_type = knot_type
        temp_val = delete1
      
        while left_side_type == knot_type:

            if delete1 + delete2 > chain_len - 3:

                return [knot_type, delete1, delete2] 

            cut = atoms[delete1:-delete2]

            left_side_type = tp.alexander(cut, closure = closure_type, max_cross = max_cross_param)
                

            if left_side_type == knot_type:
               
                delete1 += margin


            else:

                at = 1
                    
                while at <= atoms_to_check and left_side_type != knot_type:

                    cut = atoms[delete1+at*margin:-delete2]

                    if len(cut) == 0:

                        break

                    else:

                        left_side_type = tp.alexander(cut, closure = closure_type, max_cross = max_cross_param)

                        at += 1

                if left_side_type == knot_type:

                    delete1 += margin

                else:

                    delete1 -= margin
                    final_left = delete1
        delete1 = temp_val
         
        while right_side_type == knot_type:

            if delete1 + delete2 > chain_len - 3:

                return [knot_type, delete1, delete2] 
           
            cut = atoms[delete1:-delete2]
            
            right_side_type = tp.alexander(cut, closure = closure_type, max_cross = max_cross_param)
                

            if right_side_type == knot_type:
               
                delete2 += margin


            else:

                at = 1
                    
                while at <= atoms_to_check and right_side_type != knot_type:

                    cut = atoms[delete1:-delete2-at*margin]

                    if len(cut) == 0:

                        break

                    else:

                        right_side_type = tp.alexander(cut, closure = closure_type, max_cross = max_cross_param)

                        at += 1

                if right_side_type == knot_type:

                    delete2 += margin

                else:

                    delete2 -= margin
                    final_right = delete2

        if final_left < 0:
            final_left = 0
        if final_right < 1:
            final_right = 1

        return [knot_type, final_left, final_right]
    else:

        print('bad closure type')



def PrecisionCore(file, margin = 10, delete1 = 1, delete2 = 1, rollback = 2, atoms_to_check = 2, closure_type = 2, knot_cutoff = 0.42, alex_tries = 200, max_cross_param = 50):
    
    '''
    Finds the core of a knot with a cutting margin of user's choice and then runs again with a margin of 1 to ensure the most precise result.
    Takes a list of coordinates of a protein chain and returns:
    [None], if unknot or unknown
    list containing the knot type of a protein chain and it's range ex. ['3_1', 30, 20] means a trefoil knot with a knot core ranging from the 31st atom to the 20th atom from the end.
    '''

    #file - input list from file_to_list() function
    #margin - how many atoms we are cutting at a given time
    #delete1 - starting index for cutting (ex. if we don't want to start from the beginning we can put delete1 = 49 and that starts the procedure from the 50th atom)
    #delete2 - same as above except it's counted from the back (delete2 = 30 means we start from the 30th atom counting from the back)
    #rollback - when we hit the point where there are no more knots, we can go back a set amount diagonally 
    #atoms_to_check - set a number of atoms to check after finding an unknot, so that it confirms there are really no more knots after our current point
    #closure_type - sets the 'closure' parameter for the tp.alexander() function
    #knot_cutoff - set the threshold value of what we count as a knot, ex. the value returned by alexander needs to be higher than our cutoff to be recognized as a knot (only on closure_type = 2)
    #alex_tries - sets the 'tries' parameter for tp.alexander() function
    #max_cross - sets the 'max_cross' parameter for tp.alexander() function
    a = KnotCore(file, margin, delete1, delete2, rollback, atoms_to_check , closure_type , knot_cutoff, alex_tries, max_cross_param)
    if margin == 1:
        return a
    elif type(a) is not list:
        return [None]
    else:
        
        b = KnotCore(file, 1, a[1], a[2], rollback, atoms_to_check , closure_type , knot_cutoff, alex_tries, max_cross_param)

        return b


