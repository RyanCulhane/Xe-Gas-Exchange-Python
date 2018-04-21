# Import smtplib for the actual sending function
import smtplib, os
from email.MIMEMultipart import MIMEMultipart
from email.MIMEText import MIMEText
from email.MIMEBase import MIMEBase
from email import encoders

# Subject_ID = '002014'
# data_dir = '/home/ziyiw/Patients/'

def send_email(data_dir, Subject_ID):

    me = 'xemri.test@gmail.com'
    you = ["ziyi.wang@duke.edu","wangziyi1994@gmail.com"]

    msg = MIMEMultipart()
    msg['Subject'] = 'Summary for Subject '+Subject_ID #% textfile
    msg['From'] = me
    msg['To'] = ", ".join(you)

    # insert body text
    fp = open(os.path.dirname(__file__)+'/email.txt', 'rb')
    msg.attach(MIMEText(fp.read(), 'plain'))
    fp.close()


    # encoding and attach a ppt
    filename = "report_"+Subject_ID+".pptx"
    attachment = open(data_dir+'/'+filename, "rb")

    part = MIMEBase('application', 'octet-stream')
    part.set_payload((attachment).read())
    encoders.encode_base64(part)
    part.add_header('Content-Disposition', "attachment; filename= %s" % filename)

    msg.attach(part)

    # Send the message via Gmail server
    s = smtplib.SMTP(host='smtp.gmail.com', port=587)
    # initiate conversation
    s.ehlo()
    # use advanced security protocol
    s.starttls()
    # log in with credential
    s.login('xemri.test@gmail.com','4.teamxenon')
    s.sendmail(me, you, msg.as_string())
    s.quit()
