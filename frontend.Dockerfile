# build environment
FROM node:latest as builder
RUN mkdir /app
WORKDIR /app
ENV PATH /app/node_modules/.bin:$PATH
COPY gui/package.json /app/package.json
COPY gui/yarn.lock /app/yarn.lock
RUN yarn
COPY gui /app
RUN yarn build

# production environment
FROM nginx:1.13.9-alpine
COPY --from=builder /app/build /app/nomadxt
CMD ["nginx", "-g", "daemon off;"]
