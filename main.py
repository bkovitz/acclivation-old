from climber import Climber
from data import Data
import math
import sys



#returns a list of climbers of size numClimb with initial fitness
def makeClimbers(stepSize, data, numClimb):
    storage = []
    for i in range(0, len(data.table), len(data.table)//numClimb):
        (startx, starty, startFit) = data.table[i]
        climber = Climber(startx, starty, stepSize)
        climber.setStartingFitness(data)
        storage.append(climber)
    return storage
    



    
if __name__ == '__main__':
    inputFilename = sys.argv[1]
    data = Data(inputFilename)    

    numClimb = int(sys.argv[2]) # Put 157: this scatters them so they don't
                                # all show up on y=x

    climbFilename = inputFilename + '-climbers-'+ str(numClimb)
    stepSize = .005
    
    
    climbers = makeClimbers(stepSize, data, numClimb)
    with open(climbFilename, 'w') as climbFile:
        for climber in climbers:
            while (not climber.isDone()):
                climber.takeStep(data)
            #when finished, write to file
            climbFile.write(str(climber) + '\n')

    
    
            


    
 


