function R = QtoC(qe)
q1 = qe(1);q2 = qe(2);q3 = qe(3);q4 = qe(4);
R = [q1^2-q2^2-q3^2+q4^2  2*(q1*q2+q3*q4)  2*(q1*q3-q2*q4);
     2*(q1*q2-q3*q4)  q2^2-q1^2-q3^2+q4^2  2*(q3*q2+q1*q4);
     2*(q1*q3+q2*q4)  2*(q3*q2-q1*q4)  q3^2-q1^2-q2^2+q4^2];
