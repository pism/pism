const double __attribute__ unused alpha = 0.5*pow(200*pow(secpera,-1)*(-0.5*(x+0.5)*pow(pow(y,2)+1,-1)-0.6*(1-x)*pow(sec(y),2)-0.4)+200.0*pow(secpera,-1)*x*(1.4-y)*(y+1.4),2)+20000.0*pow(secpera,-2)*pow(-0.5*(x+0.5)*pow(pow(y,2)+1,-1)-0.6*(1-x)*pow(sec(y),2)-0.4,2)+0.25*pow(200*pow(secpera,-1)*(0.6*tan(y)-0.5*atan(y))-100.0*pow(secpera,-1)*(pow(x,2)+1)*(y+1.4)+100.0*pow(secpera,-1)*(pow(x,2)+1)*(1.4-y),2)+20000.0*pow(secpera,-2)*pow(x,2)*pow(1.4-y,2)*pow(y+1.4,2); 
const double __attribute__ unused alpha_x = 1.0*(200*pow(secpera,-1)*(0.6*pow(sec(y),2)-0.5*pow(pow(y,2)+1,-1))+200.0*pow(secpera,-1)*(1.4-y)*(y+1.4))*(200*pow(secpera,-1)*(-0.5*(x+0.5)*pow(pow(y,2)+1,-1)-0.6*(1-x)*pow(sec(y),2)-0.4)+200.0*pow(secpera,-1)*x*(1.4-y)*(y+1.4))+40000.0*pow(secpera,-2)*(0.6*pow(sec(y),2)-0.5*pow(pow(y,2)+1,-1))*(-0.5*(x+0.5)*pow(pow(y,2)+1,-1)-0.6*(1-x)*pow(sec(y),2)-0.4)+40000.0*pow(secpera,-2)*x*pow(1.4-y,2)*pow(y+1.4,2)+0.5*(200.0*pow(secpera,-1)*x*(1.4-y)-200.0*pow(secpera,-1)*x*(y+1.4))*(200*pow(secpera,-1)*(0.6*tan(y)-0.5*atan(y))-100.0*pow(secpera,-1)*(pow(x,2)+1)*(y+1.4)+100.0*pow(secpera,-1)*(pow(x,2)+1)*(1.4-y)); 
const double __attribute__ unused alpha_y = 1.0*(200*pow(secpera,-1)*(1.0*(x+0.5)*y*pow(pow(y,2)+1,-2)-1.2*(1-x)*tan(y)*pow(sec(y),2))-200.0*pow(secpera,-1)*x*(y+1.4)+200.0*pow(secpera,-1)*x*(1.4-y))*(200*pow(secpera,-1)*(-0.5*(x+0.5)*pow(pow(y,2)+1,-1)-0.6*(1-x)*pow(sec(y),2)-0.4)+200.0*pow(secpera,-1)*x*(1.4-y)*(y+1.4))+40000.0*pow(secpera,-2)*(1.0*(x+0.5)*y*pow(pow(y,2)+1,-2)-1.2*(1-x)*tan(y)*pow(sec(y),2))*(-0.5*(x+0.5)*pow(pow(y,2)+1,-1)-0.6*(1-x)*pow(sec(y),2)-0.4)+0.5*(200*pow(secpera,-1)*(0.6*tan(y)-0.5*atan(y))-100.0*pow(secpera,-1)*(pow(x,2)+1)*(y+1.4)+100.0*pow(secpera,-1)*(pow(x,2)+1)*(1.4-y))*(200*pow(secpera,-1)*(0.6*pow(sec(y),2)-0.5*pow(pow(y,2)+1,-1))-200.0*pow(secpera,-1)*(pow(x,2)+1))-40000.0*pow(secpera,-2)*pow(x,2)*(1.4-y)*pow(y+1.4,2)+40000.0*pow(secpera,-2)*pow(x,2)*pow(1.4-y,2)*(y+1.4); 
const double __attribute__ unused u = 100.0*pow(secpera,-1)*(pow(x,2)+1)*(1.4-y)*(y+1.4); 
const double __attribute__ unused u_x = 200.0*pow(secpera,-1)*x*(1.4-y)*(y+1.4); 
const double __attribute__ unused u_y = 100.0*pow(secpera,-1)*(pow(x,2)+1)*(1.4-y)-100.0*pow(secpera,-1)*(pow(x,2)+1)*(y+1.4); 
const double __attribute__ unused u_xx = 200.0*pow(secpera,-1)*(1.4-y)*(y+1.4);
									       
const double __attribute__ unused u_xy = 200.0*pow(secpera,-1)*x*(1.4-y)-200.0*pow(secpera,-1)*x*(y+1.4); 
const double __attribute__ unused u_yy = -200.0*pow(secpera,-1)*(pow(x,2)+1); 
const double __attribute__ unused v = 200*pow(secpera,-1)*(-0.6*(1-x)*tan(y)-0.5*(x+0.5)*atan(y)-0.4*y); 
const double __attribute__ unused v_x = 200*pow(secpera,-1)*(0.6*tan(y)-0.5*atan(y)); 
const double __attribute__ unused v_y = 200*pow(secpera,-1)*(-0.5*(x+0.5)*pow(pow(y,2)+1,-1)-0.6*(1-x)*pow(sec(y),2)-0.4); 
const double __attribute__ unused v_xx = 0; 
const double __attribute__ unused v_xy = 200*pow(secpera,-1)*(0.6*pow(sec(y),2)-0.5*pow(pow(y,2)+1,-1)); 
const double __attribute__ unused v_yy = 200*pow(secpera,-1)*(1.0*(x+0.5)*y*pow(pow(y,2)+1,-2)-1.2*(1-x)*tan(y)*pow(sec(y),2)); 
const double __attribute__ unused h = 500*(0.2*(1.2-x)*pow(y,4)+3*pow(x+4,-1));
									       
const double __attribute__ unused H = 500*(0.2*(1.2-x)*pow(y,4)+3*pow(x+4,-1))+80*x; 
const double __attribute__ unused b = -80*x; 
const double __attribute__ unused h_x = 500*(-0.2*pow(y,4)-3*pow(x+4,-2)); 
const double __attribute__ unused h_y = 400.0*(1.2-x)*pow(y,3); 
const double __attribute__ unused H_x = 500*(-0.2*pow(y,4)-3*pow(x+4,-2))+80; 
const double __attribute__ unused H_y = 400.0*(1.2-x)*pow(y,3); 
