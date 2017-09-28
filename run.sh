nohup java -Xmx16000m -classpath /usr/local/sunny/newjars/cplex.jar -Djava.library.path=/usr/local/sunny/newjars -jar exp1.jar 5 5 15 5 1 > out1.txt
sleep 10
nohup java -Xmx16000m -classpath /usr/local/sunny/newjars/cplex.jar -Djava.library.path=/usr/local/sunny/newjars -jar exp2.jar 8 8 30 8 2 > out2.txt
sleep 10
nohup java -Xmx16000m -classpath /usr/local/sunny/newjars/cplex.jar -Djava.library.path=/usr/local/sunny/newjars -jar exp3.jar 10 10 40 10 2 > out3.txt
sleep 10
nohup java -Xmx16000m -classpath /usr/local/sunny/newjars/cplex.jar -Djava.library.path=/usr/local/sunny/newjars -jar exp4.jar 14 14 50 14 2 > out4.txt
sleep 10




