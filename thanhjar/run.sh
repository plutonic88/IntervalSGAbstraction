nohup java -Xmx16000m -classpath /usr/local/sunny/newjars/cplex.jar -Djava.library.path=/usr/local/sunny/newjars -jar real.jar 10 10 5000 10 > out1.txt
sleep 10
nohup java -Xmx16000m -classpath /usr/local/sunny/newjars/cplex.jar -Djava.library.path=/usr/local/sunny/newjars -jar real.jar 10 20 8000 10 > out2.txt
sleep 10
nohup java -Xmx16000m -classpath /usr/local/sunny/newjars/cplex.jar -Djava.library.path=/usr/local/sunny/newjars -jar real.jar 10 50 15000 10 > out3.txt
sleep 10
nohup java -Xmx16000m -classpath /usr/local/sunny/newjars/cplex.jar -Djava.library.path=/usr/local/sunny/newjars -jar real.jar 20 50 18000 10 > out4.txt
sleep 10

