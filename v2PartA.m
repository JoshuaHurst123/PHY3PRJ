% This program simulates Part A of The Quantum Well Laboratory Exercise from PHY3PRJ. 
% The code plots the energy of the first n eigenvalues for a well of a width w.

title('Energy Levels in a Quantum Well')

w=3e-9 %Well width in meters.
m=9.1093837*10^-31 %Mass of electron in Kg.
h=6.62607015*10^-34 %Planck's constant in J/s.
n=[1:10]' %Quantum number
nsq=(n.*n) %Squaring matrix
% ' turns a row vector into a column vector.
% * element wise multiplication

E=(h^2)*nsq/8*m*(w^2) %From solution of Schrodinger Equation (J)
plot(n,E,'o')
hold on

% Width = 3nm. Change width from 3nm to 6nm.
w=6e-9 %Width in meters.
E=(h^2)*nsq/8*m*(w^2)
plot(n,E,'*')

% Width = 6nm. Change width from 6nm to 9nm.
w=9e-9 %Width in meters.
E=(h^2)*nsq/8*m*(w^2)
plot(n,E,'.')

xlabel('n, Quantum Number') 
ylabel('E, Energy (J)') 
legend('Wide','w=6e-9m','Narrow')