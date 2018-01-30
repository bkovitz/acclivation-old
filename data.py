class Data:

    def __init__(self, filename):

        f = open(filename, "r")
        #    data = {}
        self.table = []
        for line in f:
            
            d = line.split()
            x = float(d[0])
            y = float(d[1])
            fitness = float(d[2])
            
            self.table.append([x, y, fitness])
        
        f.close()
    
        

    def __str__(self):
        return str(self.table)

   
    #looks up the fitnes value for the given coordinate location
    def getXYfit(self, coordinates):
        x = coordinates[0]
        y = coordinates[1]
        
        fitness = float("-inf")
        for line in self.table:
            if line[0] == x:
                if line[1] == y:
                    return line[2]

        if fitness == float("-inf"):
            raise ValueError('Tried to look up fitness for coordinates not found in data table!', str(coordinates))
        return fitness
                


                
    
