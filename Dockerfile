FROM fedora

# compiler
RUN dnf -y install gcc

# multi-precision libraries
RUN dnf -y install gmp-devel mpfr-devel libmpc-devel

RUN dnf -y install python3

RUN dnf -y install python3-pip

RUN mkdir /app

COPY . /app/

WORKDIR /app

# api requirement
RUN python3 -m pip install flask

# build the binaries
RUN make build

EXPOSE 5000

ARG PORT=5000

CMD ["python3","main.py"]
