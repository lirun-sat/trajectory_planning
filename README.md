void randomlyInitial(void){
		for (int i = 0; i < _whaleCount; ++i){
			_whaleSet[i]._position[0] = 0.852354;
			_whaleSet[i]._position[1] = 0.433124;
			_whaleSet[i]._position[2] = 2.03323;
			_whaleSet[i]._position[3] = 2.9251;
			_whaleSet[i]._position[4] = 1.82466;
			_whaleSet[i]._position[5] = 1.94137;
			_whaleSet[i]._position[6] = 0.966595;
			_whaleSet[i]._fitness = _fitnessFunction(_whaleSet[i]);
		}
	}
	
	void randomlyInitial(void){
		for (int i = 0; i < _studentCount; ++i){
			_studentSet[i]._position[0] = 0.852354;
			_studentSet[i]._position[1] = 0.433124;
			_studentSet[i]._position[2] = 2.03323;
			_studentSet[i]._position[3] = 2.9251;
			_studentSet[i]._position[4] = 1.82466;
			_studentSet[i]._position[5] = 1.94137;
			_studentSet[i]._position[6] = 0.966595;
			_studentSet[i]._fitness = _fitnessFunction(_studentSet[i]);
		}
	}
