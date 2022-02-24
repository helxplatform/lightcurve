FROM python:3.9.1
ADD . /python-flask
WORKDIR /python-flask
RUN rm tmp/*
RUN pip install -r requirements.txt
RUN pip install gunicorn
USER 1001
ENV MPLCONFIGDIR=/tmp
#ENTRYPOINT [ "python3", "main.py" ]
ENTRYPOINT [ "gunicorn", "--timeout", "120", "--log-level", "debug", "--threads", "8",  "--bind", "0.0.0.0:8080", "wsgi:app" ]
