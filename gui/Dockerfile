# build environment
FROM node:9.6.1 as builder
RUN mkdir /app
WORKDIR /app
ENV PATH /app/node_modules/.bin:$PATH
COPY package.json /app/package.json
RUN npm install --silent
RUN npm install react-scripts@1.1.1 -g --silent
COPY . /app
RUN npm run build

# production environment
FROM nginx:1.13.9-alpine
COPY --from=builder /app/build /app/nomadxt
EXPOSE 80
CMD ["nginx", "-g", "daemon off;"]
