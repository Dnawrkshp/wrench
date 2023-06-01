/*
	wrench - A set of modding tools for the Ratchet & Clank PS2 games.
	Copyright (C) 2019-2023 chaoticgd

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <catch2/catch_amalgamated.hpp>

#include <core/filesystem.h>
#include <instancemgr/pvar.h>

#define CPP_TEST_PASSED -1
static s32 test_lexer(const char* src, std::vector<CppTokenType>&& expected);
static void print_token(const CppToken& token);

TEST_CASE("c++ lexer" "[instancemgr]") {
	CHECK(CPP_TEST_PASSED == test_lexer(
		"char c = \"\x42\";",
		{CPP_KEYWORD, CPP_IDENTIFIER, CPP_OPERATOR, CPP_LITERAL, CPP_OPERATOR}));
	
	CHECK(CPP_TEST_PASSED == test_lexer(
		"const char* simple_str = \"simple string\";",
		{CPP_KEYWORD, CPP_KEYWORD, CPP_OPERATOR, CPP_IDENTIFIER, CPP_OPERATOR, CPP_LITERAL, CPP_OPERATOR}));
	
	CHECK(CPP_TEST_PASSED == test_lexer(
		"const char* raw_str = R\"abc(\\\"Hello World\\\"\n(Hello\\x20World))abc\";",
		{CPP_KEYWORD, CPP_KEYWORD, CPP_OPERATOR, CPP_IDENTIFIER, CPP_OPERATOR, CPP_LITERAL, CPP_OPERATOR}));
	
	CHECK(CPP_TEST_PASSED == test_lexer(
		"struct SomeStruct {int a;}",
		{CPP_KEYWORD, CPP_IDENTIFIER, CPP_OPERATOR, CPP_KEYWORD, CPP_IDENTIFIER, CPP_OPERATOR, CPP_OPERATOR}));
}

static s32 test_lexer(const char* src, std::vector<CppTokenType>&& expected) {
	std::string str(src);
	std::vector<CppToken> tokens;
	const char* result = eat_cpp_file(tokens, &str[0]);
	if(result != CPP_NO_ERROR) {
		UNSCOPED_INFO(stringf("error: %s\n(end of error)\n", result));
	}
	for(size_t i = 0; i < std::max(tokens.size(), expected.size()); i++) {
		print_token(tokens[i]);
		if(i >= tokens.size()) return -3;
		if(i >= expected.size()) return -4;
		if(tokens[i].type != expected[i]) {
			return (s32) i;
		}
	}
	return CPP_TEST_PASSED;
}

static void print_token(const CppToken& token) {
	switch(token.type) {
		case CPP_COMMENT: {
			std::string str(token.str_begin, token.str_end);
			UNSCOPED_INFO(stringf("comment %s\n", str.c_str()));
			break;
		}
		case CPP_IDENTIFIER: {
			std::string str(token.str_begin, token.str_end);
			UNSCOPED_INFO(stringf("identifier %s\n", str.c_str()));
			break;
		}
		case CPP_KEYWORD: {
			for(s32 i = 0; i < CPP_KEYWORD_COUNT; i++) {
				if(CPP_KEYWORDS[i].keyword == token.keyword) {
					UNSCOPED_INFO(stringf("keyword %s\n", CPP_KEYWORDS[i].string));
					break;
				}
			}
			break;
		}
		case CPP_LITERAL: {
			std::string str(token.str_begin, token.str_end);
			UNSCOPED_INFO(stringf("literal %s\n", str.c_str()));
			break;
		}
		case CPP_OPERATOR: {
			for(s32 i = 0; i < CPP_OPERATOR_COUNT; i++) {
				if(CPP_OPERATORS[i].op == token.op) {
					UNSCOPED_INFO(stringf("operator %s\n", CPP_OPERATORS[i].string));
					break;
				}
			}
			break;
		}
	}
}

static bool test_parser(const char* src, PvarType&& expected);
static bool compare_pvar_types(const PvarType& lhs, const PvarType& rhs);

TEST_CASE("c++ parser" "[instancemgr]") {
	CHECK(test_parser(
		"struct SomeVars { int array_of_ints[5]; };",
		[]() {
			PvarType type(PTD_STRUCT_OR_UNION);
			type.name = "SomeVars";
			PvarType& field = type.struct_or_union.fields.emplace_back(PTD_ARRAY);
			field.name = "array_of_ints";
			field.array.element_count = 5;
			field.array.element_type = std::make_unique<PvarType>(PTD_BUILT_IN);
			field.array.element_type->built_in = PvarBuiltIn::INT;
			return type;
		}()
	));
}

static bool test_parser(const char* src, PvarType&& expected) {
	std::string str(src);
	std::vector<CppToken> tokens;
	const char* result = eat_cpp_file(tokens, &str[0]);
	if(result != CPP_NO_ERROR) UNSCOPED_INFO(stringf("error: %s\n(end of error)\n", result));
	std::vector<PvarType> types = parse_pvar_types(tokens);
	if(types.size() != 1) {
		return false;
	}
	return compare_pvar_types(types[0], expected);
}

static bool compare_pvar_types(const PvarType& lhs, const PvarType& rhs) {
	if(lhs.name != rhs.name) { UNSCOPED_INFO("name"); return false; }
	if(lhs.offset != rhs.offset) { UNSCOPED_INFO("offset"); return false; }
	if(lhs.size != rhs.size) { UNSCOPED_INFO("size"); return false; }
	if(lhs.alignment != rhs.alignment) { UNSCOPED_INFO("alignment"); return false; }
	if(lhs.descriptor != rhs.descriptor) { UNSCOPED_INFO("descriptor"); return false; }
	switch(lhs.descriptor) {
		case PTD_ARRAY: {
			if(lhs.array.element_count != rhs.array.element_count) { UNSCOPED_INFO("array.element_count"); return false; }
			REQUIRE((lhs.array.element_type.get() && rhs.array.element_type.get()));
			bool comp_result = compare_pvar_types(*lhs.array.element_type.get(), *rhs.array.element_type.get());
			if(!comp_result) { UNSCOPED_INFO("array.element_type"); return false; }
			break;
		}
		case PTD_BUILT_IN: {
			if(lhs.built_in != rhs.built_in) { UNSCOPED_INFO("built_in"); return false; }
			break;
		}
		case PTD_STRUCT_OR_UNION: {
			if(lhs.struct_or_union.is_union != rhs.struct_or_union.is_union) { UNSCOPED_INFO("struct_or_union.is_union"); return false; }
			if(lhs.struct_or_union.fields.size() != rhs.struct_or_union.fields.size()) { UNSCOPED_INFO("struct_or_union.fields.size()"); return false; }
			for(s32 i = 0; i < (s32) lhs.struct_or_union.fields.size(); i++) {
				bool comp_result = compare_pvar_types(lhs.struct_or_union.fields[i], rhs.struct_or_union.fields[i]);
				if(!comp_result) { UNSCOPED_INFO(stringf("struct_or_union.fields[%d]", i)); return false; }
			}
			break;
		}
		case PTD_TYPE_NAME: {
			UNSCOPED_INFO("unhandled type name");
			return false;
		}
		case PTD_POINTER_OR_REFERENCE: {
			if(lhs.pointer_or_reference.is_reference != rhs.pointer_or_reference.is_reference) { UNSCOPED_INFO("pointer_or_reference.is_reference"); return false; }
			REQUIRE((lhs.pointer_or_reference.value_type.get() && rhs.pointer_or_reference.value_type.get()));
			bool comp_result = compare_pvar_types(*lhs.pointer_or_reference.value_type.get(), *rhs.pointer_or_reference.value_type.get());
			if(!comp_result) { UNSCOPED_INFO("pointer_or_reference.value_type"); return false; }
			break;
		}
	}
	return true;
}
