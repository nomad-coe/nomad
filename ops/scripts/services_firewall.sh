# rabbit
firewall-cmd --permanent --add-port=5672/tcp

# mongo
firewall-cmd --permanent --add-port=27017/tcp
firewall-cmd --permanent --add-port=37017/tcp

# raw api
firewall-cmd --permanent --add-port=18001/tcp

# es
firewall-cmd --permanent --add-port=19200/tcp

# postgres
firewall-cmd --permanent --add-port=5432/tcp

# elk
firewall-cmd --permanent --add-port=29200/tcp
firewall-cmd --permanent --add-port=15601/tcp
firewall-cmd --permanent --add-port=15000/tcp

firewall-cmd --reload