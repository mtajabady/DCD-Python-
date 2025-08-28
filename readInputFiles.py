class RWfile:
    def __init__(self):
        self.hh=1
    def readFile(file_path):            # insert Sting Path of files to read , example: "D:/user/data.txt"
        fileObj = open(file_path, "r") #opens the file in read mode
        words = fileObj.read().splitlines() #puts the file into an array
        fileObj.close()
        return words

    def list_input_float(the_list , i ):    # Converting list of strings to list of float  Numbers
        disallowed_characters = "{!}"       # Remove unneeded Characters
        for character in disallowed_characters:
            the_list[i] = the_list[i].replace(character, "")

        x= []
        x_1 = list((the_list[i]).split(",")) #define the character to split the list
        for i in range(len(x_1)):
            x[len(x):] = [float(x_1[i])]
        
        return x

    def list_input_complex(the_list , i ):  # Converting list of strings to list of complex Numbers
        disallowed_characters = "{!}*"       # Remove unneeded Characters
        replace_Characters = "I"
        for character in replace_Characters:
            the_list[i] = the_list[i].replace(character, "j")
        for character in disallowed_characters:
            the_list[i] = the_list[i].replace(character, "")


        x= []
        x_1 = list((the_list[i]).split(",")) #define the character to split the list
        
        for i in range(len(x_1)):
            x[len(x):] = [complex(x_1[i])]
        
        return x

    def write_output (list_of_data , output_path):
        for i in range(len(list_of_data)):
            list_of_data[i] = str(list_of_data[i])

        with open(output_path, 'w') as f:
            for line in list_of_data:
                f.write(line)
                f.write('\n')
        return 0