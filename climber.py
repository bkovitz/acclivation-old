class Climber:
    def __init__(self, startx, starty, stepSize):
        self.startx = startx
        self.starty = starty
        self.x = startx
        self.y = starty
        self.fitness = float("-inf") #default fitness value
        self.i = 0 #i is for number of iterations completed
        self.step = stepSize
        self.xMin = -1
        self.yMin = -1
        self.xMax = 1
        self.yMax = 1
        self.done = False

    #prints the climber as: 'starting-x starting-y starting-fitness current-x current-y number-of-iterations-taken'
    def __str__(self):
        return "{0} {1} {2} {3} {4} {5} {6}".format(self.startx, self.starty, self.startFit, self.x, self.y, self.fitness, self.i)

    def isDone(self):
        return self.done
    

    def setStartingFitness(self, data):
        coordinates = [self.x, self.y]
        fitness = data.getXYfit(coordinates)
        self.startFit = fitness
        self.fitness = fitness
                    
        
    #Finds the neighbor among cardinal directions with the highest fitness for current climber
    def findBestNeighbor(self, data):

   
        north = [self.x, round(min(self.y + self.step, self.yMax), 6)]
        
        south = [self.x,  round(max(self.y - self.step, self.yMin), 6)]

        east = [round(min(self.x + self.step, self.xMax), 6), self.y]

        west = [round(max(self.x - self.step, self.xMin), 6), self.y]

        stay = [self.x, self.y]
       

        fitnesses = {}
        try:
            northFit = data.getXYfit(north)
        except ValueError as err:
            print(err.args())

        try:
            southFit = data.getXYfit(south)
        except ValueError as err:
            print(err.args())

        try:    
            eastFit = data.getXYfit(east)
        except ValueError as err:
            print(err.args())
            
        try:
            westFit = data.getXYfit(west)
        except ValueError as err:
            print(err.args())
            
        stayFit = self.fitness

        fitnesses[northFit] = north
        fitnesses[southFit] = south
        fitnesses[eastFit] = east
        fitnesses[westFit]= west
        fitnesses[stayFit] = stay


        fitlist = list(fitnesses.keys())
        fitlist.sort()

        bestfit = fitlist[-1]



        best = fitnesses[bestfit]
        best.append(bestfit)


        return best
        
    #takes a step to the best neighbor. If best option is to stay put then the climber stops
    def takeStep(self, data):

        

        best = self.findBestNeighbor(data)
        
        newX = best[0]
        newY = best[1]
        newFit = best[2]

        if (self.x == newX and self.y == newY):
            self.done = True
        else:    
            self.x = newX
            self.y = newY
            self.fitness = newFit
            self.i += 1
        
        
        

        
        
        
                
                
        




    
