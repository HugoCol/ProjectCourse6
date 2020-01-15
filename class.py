class auto:
    def __init__(self, v):
        self.__snelheid = v

    def setSnelheid(self, v):
        if v>0 and v<200:
            self.__snelheid = v
        else: print('de snelheid klopt niet')

    def getSnelheid(self):
        return self.__snelheid

