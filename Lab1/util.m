clear , close all;

global PATH
PATH = 'G:/DIP/Lab1/';


I = imread(strcat(PATH, 'uni.jpg'));   
I = rgb2gray(I);
imwrite(I, strcat(PATH, 'Uni.png'));