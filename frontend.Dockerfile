# build environment
FROM node:latest as builder
RUN mkdir -p /nomad/app
WORKDIR /nomad/app
ENV PATH /nomad/app/node_modules/.bin:$PATH
COPY gui/package.json /nomad/app/package.json
COPY gui/yarn.lock /nomad/app/yarn.lock
RUN yarn
COPY gui /nomad/app
COPY .git /nomad

RUN yarn build

# production environment
FROM nginx:1.13.9-alpine
COPY --from=builder /nomad/app/build /app/nomadxt
CMD ["nginx", "-g", "daemon off;"]
