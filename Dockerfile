
# Dockerfile 04/2021
# Python environment for exporting data analysis for processing on remote systems. 

FROM python:3.8.5
MAINTAINER Alex Perkins <a.j.p.perkins@sms.ed.ac.uk>
WORKDIR /app
COPY . /app
RUN python -m pip install -r requirements.txt
#CMD ["python", "a280.py"]