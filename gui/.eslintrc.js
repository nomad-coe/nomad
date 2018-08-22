module.exports = {
    "extends": ["standard", "plugin:react/recommended"],
    "parser": "babel-eslint",
    "parserOptions": {
        "ecmaFeatures": {
            "jsx": true,
            "modules": true
        }
    },
    "globals": {
        "fetch": false
    },
    "rules": {
        "space-before-function-paren": ["error", "never"]
    }
}