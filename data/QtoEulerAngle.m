function EA = QtoEulerAngle(Q)
%% q = [q1 q2 q3 q4]';qv = [q1 q2 q3]';q0 = q4;
EA =180/pi*[atan2(2*(Q(4)*Q(1)+Q(2)*Q(3)),(Q(4)^2-Q(1)^2-Q(2)^2+Q(3)^2));
            asin(2*Q(4)*Q(2)-2*Q(1)*Q(3));
            atan2(2*(Q(4)*Q(3)+Q(1)*Q(2)),(Q(4)^2+Q(1)^2-Q(2)^2-Q(3)^2))];