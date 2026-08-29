FROM python:3.10.5-slim-buster

WORKDIR /app

COPY ./scripts/motif_search/ .

COPY ./prody ./prody

RUN apt update -y && \
    apt install build-essential -y

RUN pip install --upgrade pip setuptools wheel

RUN pip install -r requirements.txt

RUN rm requirements.txt

CMD ["flask", "run", "--host", "0.0.0.0", "--port", "5000"]

