# build environment
FROM node:latest as builder
RUN mkdir /app
WORKDIR /app
ENV PATH /app/node_modules/.bin:$PATH
COPY package.json /app/package.json
COPY yarn.lock /app/yarn.lock
RUN yarn
COPY . /app
RUN yarn build

# production environment
FROM nginx:1.13.9-alpine
COPY --from=builder /app/build /app/nomadxt
CMD ["nginx", "-g", "daemon off;"]
